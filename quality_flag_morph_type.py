# -*- coding: utf-8 -*-
"""
Look at how the quality flag filters down our samples

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

from vars import phangs_folder, output_folder, pattern_speed_version, alma_galaxies, galaxy_table, plot_folder

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

corot_table = Table.read(output_folder + 'pattern_speed_table_' + pattern_speed_version + '.fits')
# print(corot_table)

total_alma_morph = []
quality_alma_morph = []

total_muse_mass_morph = []
quality_muse_mass_morph = []

total_muse_ha_morph = []
quality_muse_ha_morph = []

x = np.arange(-8, 11)
x_min = np.min(x)
x_max = np.max(x)

for galaxy in alma_galaxies:

    phangs_row = galaxy_table[galaxy_table['name'] == galaxy]
    morph = phangs_row['morph_t'][0]

    corot_row = corot_table[corot_table['GALAXY'] == galaxy.upper()]

    try:
        alma_q = corot_row['OM_P_ALMA_QUAL'][0]
    except IndexError:
        continue

    total_alma_morph.append(morph)
    if alma_q in [1, 2]:
        quality_alma_morph.append(morph)

for galaxy in alma_galaxies:

    phangs_row = galaxy_table[galaxy_table['name'] == galaxy]
    morph = phangs_row['morph_t'][0]

    corot_row = corot_table[corot_table['GALAXY'] == galaxy.upper()]

    try:
        muse_q = corot_row['OM_P_MUSE_MASS_QUAL'][0]
    except IndexError:
        continue

    if not np.isnan(muse_q):
        total_muse_mass_morph.append(morph)
        if muse_q in [1, 2]:
            quality_muse_mass_morph.append(morph)

for galaxy in alma_galaxies:

    phangs_row = galaxy_table[galaxy_table['name'] == galaxy]
    morph = phangs_row['morph_t'][0]

    corot_row = corot_table[corot_table['GALAXY'] == galaxy.upper()]

    try:
        muse_q = corot_row['OM_P_MUSE_HA_QUAL'][0]
    except IndexError:
        continue

    if not np.isnan(muse_q):
        total_muse_ha_morph.append(morph)
        if muse_q in [1, 2]:
            quality_muse_ha_morph.append(morph)

plt.figure(figsize=(4.5, 5))

plt.subplot(3, 1, 3)
plt.hist(total_alma_morph, bins=x, color='grey', label=r'All CO: $N=%s$' % len(total_alma_morph))
plt.hist(quality_alma_morph, bins=x, color='b', label=r'$Q = 1, 2$: $N=%s$' % (len(quality_alma_morph)))

# plt.xticks([])
plt.xlim([x_min, x_max])

plt.xlabel(r'T$_\mathrm{Hubble}$')
plt.ylabel(r'$N_\mathrm{gal}$')

plt.legend(loc='upper left', frameon=False)

plt.subplot(3, 1, 1)
plt.hist(total_muse_mass_morph, bins=x, color='grey', label=r'All $M_\ast$: $N=%s$' % len(total_muse_mass_morph))
plt.hist(quality_muse_mass_morph, bins=x, color='r', label='$Q = 1, 2$: $N=%s$' % len(quality_muse_mass_morph))

plt.xticks([])
plt.xlim([x_min, x_max])

plt.ylabel(r'$N_\mathrm{gal}$')

plt.legend(loc='upper left', frameon=False)

plt.subplot(3, 1, 2)
plt.hist(total_muse_ha_morph, bins=x, color='grey', label=r'All H$\alpha$: $N=%s$' % len(total_muse_ha_morph))
plt.hist(quality_muse_ha_morph, bins=x, color='cyan', label='$Q = 1, 2$: $N=%s$' % len(quality_muse_ha_morph))

plt.legend(loc='upper left', frameon=False)

plt.xticks([])
plt.xlim([x_min, x_max])

plt.ylabel(r'$N_\mathrm{gal}$')

plt.subplots_adjust(wspace=0, hspace=0)

# plt.show()

plt.savefig(plot_folder + 'quality_flag_morph_type.png',
            bbox_inches='tight')
plt.savefig(plot_folder + 'quality_flag_morph_type.pdf',
            bbox_inches='tight')

print('Complete!')
