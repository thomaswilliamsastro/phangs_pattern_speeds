# -*- coding: utf-8 -*-
"""
Published literature pattern speeds comparison to this work

@author: Tom Williams
"""

import os

import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from scipy.stats import median_absolute_deviation

from vars import phangs_folder, output_folder, pattern_speed_version, plot_folder

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

# Dictionary here is galaxy name, pattern speed, error, distance (Mpc), inclination (deg). Comment is arXiv number or
# ADS bibcode.

rad_conv = np.radians(1 / 3600)

pattern_speeds = {
    'NGC0253': [61.3, 0, 3.47, 74, '2014A%26A...567A..86I'],
    'NGC0628': [41.8, 1, 9.93, 7, '2014ApJ...790..118M'],
    'NGC1097': [36, 2, 14.5, 35, '2014MNRAS.438..971P'],
    'NGC1300': [20, 0, 20, 35, '1996A%26A...313..733L'],
    'NGC1365': [2.4 / (18.2e3 * rad_conv), 0.5 / (18.2e3 * rad_conv), 18.2, 41, '2016ApJ...826....2S'],
    'NGC1433': [15.8, 0, 11.6, 33, '2008AJ....136..300T'],
    'NGC1512': [50, 0, 9.5, 51, '2009MNRAS.400.1749K'],
    'NGC1566': [33, 0, 17.4, 26, '2005astro.ph..9708K'],
    'NGC1672': [30, 0, 15.1, 40, '1999ApJ...512..623D'],
    'NGC2903': [40, 0, 7.3, 61, '2009PASJ...61..441H'],
    'NGC3627': [50, 0, 11.1, 65, '2004ApJ...614..142R'],
    'NGC4254': [26, 10, 16.1, 34, '2004PASJ...56L..45E'],
    'NGC4303': [46, 6, 16.1, 25, '2002ApJ...575..826S'],
    'NGC4321': [20, 1, 16.1, 31.7, '2005ApJ...632..253H'],
    'NGC4579': [52.7, 6.6, 16.5, 38.7, '2019PASJ...71S..16S'],
    'NGC4596': [52, 13, 15.7, 38, '1999MNRAS.306..926G'],
    'NGC5248': [13, 0, 15.4, 50, '2001MNRAS.324..891L'],
    'NGC6300': [27, 8, 14.3, 52, '1996ApJ...460..665R'],
}

corot_table = Table.read(output_folder + 'pattern_speed_table_' + pattern_speed_version + '.fits')
# print(corot_table.colnames)

for row in corot_table:

    if not np.isnan(row['OM_P_ALMA']) and not np.isnan(row['OM_P_MUSE_HA']):

        if row['GALAXY'] not in pattern_speeds.keys():

            pattern_speeds[row['GALAXY']] = [np.nan, np.nan, np.nan, np.nan, '']

max_len = len(max([pattern_speeds[key][4] for key in pattern_speeds.keys()], key=len))

# Add in columns to the tables for the literature speed errors
try:
    corot_table.replace_column(col=np.nan * np.zeros(len(corot_table)), name='LITERATURE_OM_P')
    corot_table.replace_column(col=np.nan * np.zeros(len(corot_table)), name='LITERATURE_OM_P_ERR')
    corot_table.replace_column(col=[' ' * max_len] * len(corot_table), name='LITERATURE_OM_P_REF')
except ValueError:
    corot_table.add_column(col=np.nan * np.zeros(len(corot_table)), name='LITERATURE_OM_P')
    corot_table.add_column(col=np.nan * np.zeros(len(corot_table)), name='LITERATURE_OM_P_ERR')
    corot_table.add_column(col=[' ' * max_len] * len(corot_table), name='LITERATURE_OM_P_REF')

corot_table['LITERATURE_OM_P'].unit = u.km / u.s / u.kpc
corot_table['LITERATURE_OM_P_ERR'].unit = u.km / u.s / u.kpc

# See how these compare to the PHANGS values

lit_pattern_speeds = []
lit_pattern_speeds_err = []
lit_pattern_speeds_y = []

lit_ism_pattern_speeds = []
lit_ism_pattern_speeds_y = []

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

galaxies_mad = []
galaxies = []

pos = 0

for i, key in enumerate(pattern_speeds.keys()):
    # Pull out the working distances/inclinations from the corotation table

    phangs_idx = np.where(corot_table['GALAXY'] == key)
    phangs_dist = corot_table[phangs_idx]['DIST'][0]
    phangs_incl = corot_table[phangs_idx]['INCL'][0]

    galaxy_mad = []

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
    lit_ref = pattern_speeds[key][4]

    # Account for differences in inclination

    lit_pattern_speed /= np.sin(np.radians(phangs_incl)) / np.sin(np.radians(lit_incl))
    lit_pattern_speed_err /= np.sin(np.radians(phangs_incl)) / np.sin(np.radians(lit_incl))

    # Account for differences in distance

    lit_pattern_speed /= phangs_dist / lit_dist
    lit_pattern_speed_err /= phangs_dist / lit_dist

    corot_table_idx = np.where(corot_table['GALAXY'] == key)[0][0]
    corot_table['LITERATURE_OM_P'][corot_table_idx] = lit_pattern_speed
    if lit_pattern_speed_err != 0:
        corot_table['LITERATURE_OM_P_ERR'][corot_table_idx] = lit_pattern_speed_err
    corot_table['LITERATURE_OM_P_REF'][corot_table_idx] = lit_ref

    if not np.isnan(muse_mass_pattern_speed_q) and not np.isnan(lit_pattern_speed):

        galaxies.append(key)
        print(key, muse_mass_pattern_speed, lit_pattern_speed)

        lit_pattern_speeds.append(lit_pattern_speed)
        lit_pattern_speeds_err.append(lit_pattern_speed_err)
        lit_pattern_speeds_y.append(pos + 1.3)

        # Highlight the problematic ones

        if key in ['NGC1365', 'NGC3627', 'NGC4321']:

            lit_ism_pattern_speeds.append(lit_pattern_speed)
            lit_ism_pattern_speeds_y.append(pos + 1.3)

        if alma_pattern_speed_q in [1, 2]:
            alma_pattern_speeds.append(alma_pattern_speed)
            alma_pattern_speeds_err_up.append(alma_pattern_speed_err_up)
            alma_pattern_speeds_err_down.append(alma_pattern_speed_err_down)
            alma_pattern_speeds_y.append(pos + 0.7)
            # galaxy_mad.append(alma_pattern_speed)
        else:
            alma_pattern_speeds_bad.append(alma_pattern_speed)
            alma_pattern_speeds_err_up_bad.append(alma_pattern_speed_err_up)
            alma_pattern_speeds_err_down_bad.append(alma_pattern_speed_err_down)
            alma_pattern_speeds_y_bad.append(pos + 0.7)

        if muse_mass_pattern_speed_q in [1, 2]:
            muse_mass_pattern_speeds.append(muse_mass_pattern_speed)
            muse_mass_pattern_speeds_err_up.append(muse_mass_pattern_speed_err_up)
            muse_mass_pattern_speeds_err_down.append(muse_mass_pattern_speed_err_down)
            muse_mass_pattern_speeds_y.append(pos + 1)
            galaxy_mad.append(muse_mass_pattern_speed)
        else:
            muse_mass_pattern_speeds_bad.append(muse_mass_pattern_speed)
            muse_mass_pattern_speeds_err_up_bad.append(muse_mass_pattern_speed_err_up)
            muse_mass_pattern_speeds_err_down_bad.append(muse_mass_pattern_speed_err_down)
            muse_mass_pattern_speeds_y_bad.append(pos + 1)

        if muse_ha_pattern_speed_q in [1, 2]:
            muse_ha_pattern_speeds.append(muse_ha_pattern_speed)
            muse_ha_pattern_speeds_err_up.append(muse_ha_pattern_speed_err_up)
            muse_ha_pattern_speeds_err_down.append(muse_ha_pattern_speed_err_down)
            muse_ha_pattern_speeds_y.append(pos + 1.3)
            # galaxy_mad.append(muse_ha_pattern_speed)
        else:
            muse_ha_pattern_speeds_bad.append(muse_ha_pattern_speed)
            muse_ha_pattern_speeds_err_up_bad.append(muse_ha_pattern_speed_err_up)
            muse_ha_pattern_speeds_err_down_bad.append(muse_ha_pattern_speed_err_down)
            muse_ha_pattern_speeds_y_bad.append(pos + 1.3)

        if len(galaxy_mad) > 0 and not np.isnan(lit_pattern_speed) and key != 'NGC3627':
            galaxy_mad = np.array(galaxy_mad)
            galaxy_mad -= lit_pattern_speed
            galaxy_mad /= lit_pattern_speed

            galaxies_mad.extend(galaxy_mad)

        pos += 1

corot_table.write(output_folder + 'pattern_speed_table_' + pattern_speed_version + '.fits', overwrite=True)

# Calculate the MAD for the whole sample

print(galaxies_mad)
# no

galaxies_mad = np.array(galaxies_mad)
print('Percentage MAD: %.2f' % (median_absolute_deviation(galaxies_mad)*100))

plt.figure(figsize=(4, 0.25 * len(pattern_speeds.keys())))

plt.errorbar(lit_pattern_speeds, lit_pattern_speeds_y,
             xerr=lit_pattern_speeds_err, c='gray', ls='none', marker='o',
             label='Literature')
plt.scatter(lit_ism_pattern_speeds, lit_ism_pattern_speeds_y, c='none', edgecolor='lime', lw=2, zorder=99)

plt.errorbar(muse_mass_pattern_speeds, muse_mass_pattern_speeds_y,
             xerr=[muse_mass_pattern_speeds_err_down, muse_mass_pattern_speeds_err_up],
             c='r', ls='none', marker='o', label=r'MUSE M$_\ast$ (This work)')
plt.errorbar(muse_mass_pattern_speeds_bad, muse_mass_pattern_speeds_y_bad,
             xerr=[muse_mass_pattern_speeds_err_down_bad, muse_mass_pattern_speeds_err_up_bad],
             c='r', ls='none', marker='o', fillstyle='none')

plt.errorbar(alma_pattern_speeds, alma_pattern_speeds_y,
             xerr=[alma_pattern_speeds_err_down, alma_pattern_speeds_err_up],
             c='b', ls='none', marker='o', label='ALMA CO (This work)')
plt.errorbar(alma_pattern_speeds_bad, alma_pattern_speeds_y_bad,
             xerr=[alma_pattern_speeds_err_down_bad, alma_pattern_speeds_err_up_bad],
             c='b', ls='none', marker='o', fillstyle='none')

# plt.errorbar(muse_ha_pattern_speeds, muse_ha_pattern_speeds_y,
#              xerr=[muse_ha_pattern_speeds_err_down, muse_ha_pattern_speeds_err_up],
#              c='cyan', ls='none', marker='o', label=r'MUSE H$\alpha$ (This work)')
# plt.errorbar(muse_ha_pattern_speeds_bad, muse_ha_pattern_speeds_y_bad,
#              xerr=[muse_ha_pattern_speeds_err_down_bad, muse_ha_pattern_speeds_err_up_bad],
#              c='cyan', ls='none', marker='o', fillstyle='none')

for i in np.arange(len(galaxies)):
    plt.axhline(i + 1 - 0.5, c='k', ls='--')
    plt.axhline(i + 1 + 0.5, c='k', ls='--')

plt.yticks(np.arange(len(galaxies))+1, galaxies)

plt.xlim([0, 100])
plt.ylim([0, np.max(lit_pattern_speeds_y) + 3])

plt.xlabel(r'$\Omega_\mathrm{P}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')

plt.legend(loc='upper right', frameon=False)

plt.tight_layout()

# plt.show()

plt.savefig(plot_folder + 'literature_comparison.png',
            bbox_inches='tight')
plt.savefig(plot_folder + 'literature_comparison.pdf',
            bbox_inches='tight')

print('Complete!')
