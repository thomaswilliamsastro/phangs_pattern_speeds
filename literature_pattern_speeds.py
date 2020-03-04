# -*- coding: utf-8 -*-
"""
Published literature pattern speeds comparison to this work

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

from vars import phangs_folder, output_folder, corot_version, plot_folder

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

# Dictionary here is galaxy name, pattern speed, error, distance (Mpc), inclination (deg). Comment is arXiv number or
# ADS bibcode.

rad_conv = np.radians(1 / 3600)

pattern_speeds = {
    'NGC0253': [61.3, 0, 3.47, 74],  # 1405.7301
    'NGC0628': [41.8, 1, 9.93, 7],  # 1406.4561
    'NGC1097': [36, 2, 14.5, 35],  # 1311.2953
    'NGC1300': [20, 0, 20, 35],  # 1996A%26A...313..733L
    'NGC1365': [2.4 / (18.2e3 * rad_conv), 0.5 / (18.2e3 * rad_conv), 18.2, 41],  # 1606.04572, converted to km/s/kpc
    'NGC1433': [15.8, 0, 11.6, 33],  # 0804.3356, converted to km/s/kpc
    'NGC1512': [50, 0, 9.5, 51],  # 0908.4128
    'NGC1566': [33, 0, 17.4, 26],  # 0509708
    'NGC1672': [30, 0, 15.1, 40],  # 1999ApJ...512..623D
    'NGC2903': [40, 0, 7.3, 61],  # 2009PASJ...61..441H
    'NGC3627': [50, 0, 11.1, 65],  # 2004ApJ...614..142R
    'NGC4254': [26, 10, 16.1, 34],  # 0410469
    'NGC4303': [46, 6, 16.1, 25],  # 0204133
    'NGC4321': [20, 1, 16.1, 31.7],  # 0505384
    'NGC4579': [52.7, 6.6, 16.5, 38.7],  # 1901.00640
    'NGC4596': [52, 13, 15.7, 38],  # 1999MNRAS.306..926G
    'NGC5248': [13, 0, 15.4, 50],  # 0101285
    'NGC6300': [27, 8, 14.3, 52],  # 9512043
}

corot_table = Table.read(output_folder + 'pattern_speed_table_' + corot_version + '.fits')
# print(corot_table.colnames)

# See how these compare to the PHANGS values

lit_pattern_speeds = []
lit_pattern_speeds_err = []
lit_pattern_speeds_y = []

alma_pattern_speeds = []
alma_pattern_speeds_err_up = []
alma_pattern_speeds_err_down = []
alma_pattern_speeds_y = []

alma_pattern_speeds_bad = []
alma_pattern_speeds_err_up_bad = []
alma_pattern_speeds_err_down_bad = []
alma_pattern_speeds_y_bad = []

muse_mass_pattern_speeds = []
muse_mass_pattern_speeds_err_up = []
muse_mass_pattern_speeds_err_down = []
muse_mass_pattern_speeds_y = []

muse_mass_pattern_speeds_bad = []
muse_mass_pattern_speeds_err_up_bad = []
muse_mass_pattern_speeds_err_down_bad = []
muse_mass_pattern_speeds_y_bad = []

muse_ha_pattern_speeds = []
muse_ha_pattern_speeds_err_up = []
muse_ha_pattern_speeds_err_down = []
muse_ha_pattern_speeds_y = []

muse_ha_pattern_speeds_bad = []
muse_ha_pattern_speeds_err_up_bad = []
muse_ha_pattern_speeds_err_down_bad = []
muse_ha_pattern_speeds_y_bad = []

for i, key in enumerate(pattern_speeds.keys()):
    # Pull out the working distances/inclinations from the corotation table

    phangs_idx = np.where(corot_table['GALAXY'] == key)
    phangs_dist = corot_table[phangs_idx]['DIST'][0]
    phangs_incl = corot_table[phangs_idx]['INCL'][0]

    # Pull out the ALMA and MUSE pattern speeds and flags

    alma_pattern_speed = corot_table[phangs_idx]['OM_P_ALMA'][0]
    alma_pattern_speed_err_up = corot_table[phangs_idx]['OM_P_ALMA_ERR_UP'][0]
    alma_pattern_speed_err_down = corot_table[phangs_idx]['OM_P_ALMA_ERR_DOWN'][0]
    alma_pattern_speed_q = corot_table[phangs_idx]['OM_P_ALMA_QUAL'][0]

    muse_mass_pattern_speed = corot_table[phangs_idx]['OM_P_MUSE_MASS'][0]
    muse_mass_pattern_speed_err_up = corot_table[phangs_idx]['OM_P_MUSE_MASS_ERR_UP'][0]
    muse_mass_pattern_speed_err_down = corot_table[phangs_idx]['OM_P_MUSE_MASS_ERR_DOWN'][0]
    muse_mass_pattern_speed_q = corot_table[phangs_idx]['OM_P_MUSE_MASS_QUAL'][0]

    muse_ha_pattern_speed = corot_table[phangs_idx]['OM_P_MUSE_HA'][0]
    muse_ha_pattern_speed_err_up = corot_table[phangs_idx]['OM_P_MUSE_HA_ERR_UP'][0]
    muse_ha_pattern_speed_err_down = corot_table[phangs_idx]['OM_P_MUSE_HA_ERR_DOWN'][0]
    muse_ha_pattern_speed_q = corot_table[phangs_idx]['OM_P_MUSE_HA_QUAL'][0]

    # Add the literature pattern speed and associated y-position

    lit_pattern_speed = pattern_speeds[key][0]
    lit_pattern_speed_err = pattern_speeds[key][1]
    lit_dist = pattern_speeds[key][2]
    lit_incl = pattern_speeds[key][3]

    # Account for differences in inclination

    lit_pattern_speed /= np.sin(np.radians(phangs_incl)) / np.sin(np.radians(lit_incl))
    lit_pattern_speed_err /= np.sin(np.radians(phangs_incl)) / np.sin(np.radians(lit_incl))

    # Account for differences in distance

    lit_pattern_speed /= phangs_dist / lit_dist
    lit_pattern_speed_err /= phangs_dist / lit_dist

    lit_pattern_speeds.append(lit_pattern_speed)
    lit_pattern_speeds_err.append(lit_pattern_speed_err)
    lit_pattern_speeds_y.append(i + 0.9)

    if alma_pattern_speed_q in [1, 2]:
        alma_pattern_speeds.append(alma_pattern_speed)
        alma_pattern_speeds_err_up.append(alma_pattern_speed_err_up)
        alma_pattern_speeds_err_down.append(alma_pattern_speed_err_down)
        alma_pattern_speeds_y.append(i + 0.7)
    else:
        alma_pattern_speeds_bad.append(alma_pattern_speed)
        alma_pattern_speeds_err_up_bad.append(alma_pattern_speed_err_up)
        alma_pattern_speeds_err_down_bad.append(alma_pattern_speed_err_down)
        alma_pattern_speeds_y_bad.append(i + 0.7)

    if muse_mass_pattern_speed_q in [1, 2]:
        muse_mass_pattern_speeds.append(muse_mass_pattern_speed)
        muse_mass_pattern_speeds_err_up.append(muse_mass_pattern_speed_err_up)
        muse_mass_pattern_speeds_err_down.append(muse_mass_pattern_speed_err_down)
        muse_mass_pattern_speeds_y.append(i + 1.1)
    else:
        muse_mass_pattern_speeds_bad.append(muse_mass_pattern_speed)
        muse_mass_pattern_speeds_err_up_bad.append(muse_mass_pattern_speed_err_up)
        muse_mass_pattern_speeds_err_down_bad.append(muse_mass_pattern_speed_err_down)
        muse_mass_pattern_speeds_y_bad.append(i + 1.1)

    if muse_ha_pattern_speed_q in [1, 2]:
        muse_ha_pattern_speeds.append(muse_ha_pattern_speed)
        muse_ha_pattern_speeds_err_up.append(muse_ha_pattern_speed_err_up)
        muse_ha_pattern_speeds_err_down.append(muse_ha_pattern_speed_err_down)
        muse_ha_pattern_speeds_y.append(i + 1.3)
    else:
        muse_ha_pattern_speeds_bad.append(muse_ha_pattern_speed)
        muse_ha_pattern_speeds_err_up_bad.append(muse_ha_pattern_speed_err_up)
        muse_ha_pattern_speeds_err_down_bad.append(muse_ha_pattern_speed_err_down)
        muse_ha_pattern_speeds_y_bad.append(i + 1.3)

plt.figure(figsize=(4, 0.5 * len(pattern_speeds.keys()) + 1))

plt.errorbar(lit_pattern_speeds, lit_pattern_speeds_y,
             xerr=lit_pattern_speeds_err, c='gray', ls='none', marker='o',
             label='Literature')

plt.errorbar(alma_pattern_speeds, alma_pattern_speeds_y,
             xerr=[alma_pattern_speeds_err_down, alma_pattern_speeds_err_up],
             c='b', ls='none', marker='o', label='ALMA (This work)')
plt.errorbar(alma_pattern_speeds_bad, alma_pattern_speeds_y_bad,
             xerr=[alma_pattern_speeds_err_down_bad, alma_pattern_speeds_err_up_bad],
             c='b', ls='none', marker='o', fillstyle='none')

plt.errorbar(muse_mass_pattern_speeds, muse_mass_pattern_speeds_y,
             xerr=[muse_mass_pattern_speeds_err_down, muse_mass_pattern_speeds_err_up],
             c='r', ls='none', marker='o', label=r'MUSE M$_\ast$ (This work)')
plt.errorbar(muse_mass_pattern_speeds_bad, muse_mass_pattern_speeds_y_bad,
             xerr=[muse_mass_pattern_speeds_err_down_bad, muse_mass_pattern_speeds_err_up_bad],
             c='r', ls='none', marker='o', fillstyle='none')

plt.errorbar(muse_ha_pattern_speeds, muse_ha_pattern_speeds_y,
             xerr=[muse_ha_pattern_speeds_err_down, muse_ha_pattern_speeds_err_up],
             c='cyan', ls='none', marker='o', label=r'MUSE H$\alpha$ (This work)')
plt.errorbar(muse_ha_pattern_speeds_bad, muse_ha_pattern_speeds_y_bad,
             xerr=[muse_ha_pattern_speeds_err_down_bad, muse_ha_pattern_speeds_err_up_bad],
             c='cyan', ls='none', marker='o', fillstyle='none')

for i in np.arange(len(pattern_speeds.keys())):
    plt.axhline(i + 1 - 0.5, c='k', ls='--')
    plt.axhline(i + 1 + 0.5, c='k', ls='--')

plt.yticks(np.arange(len(pattern_speeds.keys()))+1, list(pattern_speeds.keys()))

plt.xlim([0, 100])
plt.ylim([0, np.max(lit_pattern_speeds_y) + 3.5])

plt.xlabel(r'$\Omega_{p, \mathrm{TW}}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')

plt.legend(loc='upper right', frameon=False)

plt.tight_layout()

# plt.show()

plt.savefig(plot_folder + 'literature_comparison.png',
            bbox_inches='tight')
plt.savefig(plot_folder + 'literature_comparison.pdf',
            bbox_inches='tight')

print('Complete!')
