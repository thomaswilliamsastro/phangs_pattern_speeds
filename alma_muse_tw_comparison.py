# -*- coding: utf-8 -*-
"""
Compare MUSE/ALMA TW pattern speeds

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

from vars import phangs_folder, alma_version, muse_version, muse_galaxies, plot_folder, output_folder, corot_version

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

# TODO: Include the Halpha

os.chdir(phangs_folder)

corot_table = Table.read(output_folder + 'pattern_speed_table_' + corot_version + '.fits')

plot_filename = plot_folder + 'muse_' + muse_version + '_'
plot_filename += 'alma_' + alma_version
plot_filename += '_pattern_speeds_comparison'

plt.figure(figsize=(4, 0.5 * len(muse_galaxies) + 1))
ax1 = plt.subplot(1, 1, 1)

frame1 = plt.gca()

position = 0

# Plot on some dummy errorbars for labels

plt.errorbar(-100, 100, xerr=3,
             fmt='o', c='r', label=r'MUSE $M_\ast$')
plt.errorbar(-100, 100, xerr=3,
             fmt='o', c='cyan', label=r'MUSE H$\alpha$')
plt.errorbar(-100, 100, xerr=3,
             fmt='o', c='b', label='ALMA')

galaxies = []

for galaxy in corot_table['GALAXY']:

    row = corot_table[corot_table['GALAXY'] == galaxy]

    muse_om_mass = row['OM_P_MUSE_MASS'][0]
    muse_om_mass_err_up = row['OM_P_MUSE_MASS_ERR_UP'][0]
    muse_om_mass_err_down = row['OM_P_MUSE_MASS_ERR_DOWN'][0]
    muse_om_mass_q = row['OM_P_MUSE_MASS_QUAL'][0]

    if np.isnan(muse_om_mass):
        continue

    muse_om_ha = row['OM_P_MUSE_HA'][0]
    muse_om_ha_err_up = row['OM_P_MUSE_HA_ERR_UP'][0]
    muse_om_ha_err_down = row['OM_P_MUSE_HA_ERR_DOWN'][0]
    muse_om_ha_q = row['OM_P_MUSE_HA_QUAL'][0]

    alma_om = row['OM_P_ALMA'][0]
    alma_om_err_up = row['OM_P_ALMA_ERR_UP'][0]
    alma_om_err_down = row['OM_P_ALMA_ERR_DOWN'][0]
    alma_om_q = row['OM_P_ALMA_QUAL'][0]

    galaxies.append(galaxy)

    if muse_om_mass_q in [1, 2]:
        fill_style = 'full'
    else:
        fill_style = 'none'

    plt.errorbar(muse_om_mass, position - 0.125,
                 xerr=np.array([[muse_om_mass_err_down, muse_om_mass_err_up]]).T,
                 fmt='o', c='r', fillstyle=fill_style)

    if muse_om_ha_q in [1, 2]:
        fill_style = 'full'
    else:
        fill_style = 'none'

    plt.errorbar(muse_om_ha, position,
                 xerr=np.array([[muse_om_ha_err_down, muse_om_ha_err_up]]).T,
                 fmt='o', c='cyan', fillstyle=fill_style)

    if alma_om_q in [1, 2]:
        fill_style = 'full'
    else:
        fill_style = 'none'

    plt.errorbar(alma_om, position + 0.125,
                 xerr=np.array([[alma_om_err_down, alma_om_err_up]]).T,
                 fmt='o', c='b', fillstyle=fill_style)

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

plt.ylim([-0.25, 0.5 * len(galaxies) + 0.75])
plt.xlim([0, 100])

ax1.yaxis.set_ticks_position('none')

plt.tight_layout()

# plt.show()

plt.savefig(plot_filename + '.png',
            bbox_inches='tight')
plt.savefig(plot_filename + '.pdf',
            bbox_inches='tight')

plt.close()

print('Complete!')
