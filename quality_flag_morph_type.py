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

from vars import phangs_folder, output_folder, corot_version, alma_galaxies, galaxy_table, plot_folder

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

corot_table = Table.read(output_folder + 'pattern_speed_table_' + corot_version + '.fits')
# print(corot_table)

total_alma_morph = []
quality_alma_morph = []

total_muse_morph = []
quality_muse_morph = []

x = np.arange(-8, 11)

for galaxy in alma_galaxies:

    phangs_row = galaxy_table[galaxy_table['NAME'] == galaxy]
    morph = phangs_row['MORPH_T'][0]

    galaxy = galaxy.strip()

    corot_row = corot_table[corot_table['GALAXY'] == galaxy]

    try:
        alma_q = corot_row['OM_P_ALMA_QUAL'][0]
    except IndexError:
        continue

    total_alma_morph.append(morph)
    if alma_q in [1, 2]:
        quality_alma_morph.append(morph)

for galaxy in alma_galaxies:

    phangs_row = galaxy_table[galaxy_table['NAME'] == galaxy]
    morph = phangs_row['MORPH_T'][0]

    galaxy = galaxy.strip()

    corot_row = corot_table[corot_table['GALAXY'] == galaxy]

    try:
        muse_q = corot_row['OM_P_MUSE_QUAL'][0]
    except IndexError:
        continue

    if muse_q != -1:
        total_muse_morph.append(morph)
        if muse_q in [1, 2]:
            quality_muse_morph.append(morph)

plt.figure(figsize=(8, 6))

plt.subplot(2, 1, 1)
plt.hist(total_alma_morph, bins=x, color='grey', label='All ALMA Pattern Speeds')
plt.hist(quality_alma_morph, bins=x, color='b', label='Quality 1 or 2')

plt.xticks([])

plt.ylabel(r'$N_\mathrm{gal}$')

plt.legend(loc='upper left', frameon=False)

plt.subplot(2, 1, 2)
plt.hist(total_muse_morph, bins=x, color='grey', label='All MUSE Pattern Speeds')
plt.hist(quality_muse_morph, bins=x, color='r', label='Quality 1 or 2')

plt.legend(loc='upper left', frameon=False)

plt.xlabel(r'T$_\mathrm{morph}$')
plt.ylabel(r'$N_\mathrm{gal}$')

plt.subplots_adjust(wspace=0, hspace=0)

plt.savefig(plot_folder + 'quality_flag_morph_type.png',
            bbox_inches='tight')
plt.savefig(plot_folder + 'quality_flag_morph_type.pdf',
            bbox_inches='tight')

print('Complete!')
