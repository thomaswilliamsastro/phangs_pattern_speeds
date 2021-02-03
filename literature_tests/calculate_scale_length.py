# -*- coding: utf-8 -*-
"""
For each galaxy, pull disc scale length

@author: Tom Williams
"""

import os

import numpy as np

from vars import phangs_folder, alma_galaxies

os.chdir(phangs_folder)

f = open('lit_tables/s4g_salo_tab_7.txt', 'r')

for i in range(44):
    f.readline()

s4g_dict = {}

s4g_galaxies = []
s4g_hr = []
s4g_bulge_total_ratio = []
s4g_flux_frac = []

for galaxy_id in range(1, 2353):
    bad_comps = 0

    line = f.readline()

    # Pull out galaxy name, number of components and quality flag
    galaxy = line[6:16].strip()
    n_comp = int(line[26:33].split('=')[-1])
    quality = int(line[35:43].split('=')[-1])

    s4g_dict[galaxy] = {'n_comp': n_comp,
                        'quality': quality}

    s4g_dict[galaxy]['flux_fracs'] = []

    for i in range(n_comp):

        line = f.readline()

        if quality not in [4, 5]:
            continue

        if line[6:16].strip() in ['Z', 'D', 'B']:

            # Pull out the bulge flux fraction

            if line[6:16].strip() == 'B':
                s4g_dict[galaxy]['b_t_ratio'] = float(line[26:33])
                continue

            # Pull out the scale length

            flux_frac = float(line[26:33])

            if line[17:25].strip() == 'expdisk':

                try:
                    if flux_frac > s4g_dict[galaxy]['flux_fracs'][-1]:
                        s4g_dict[galaxy]['h_r'] = float(line[61:67].strip())
                        s4g_dict[galaxy]['flux_fracs'].append(flux_frac)
                except IndexError:
                    s4g_dict[galaxy]['h_r'] = float(line[61:67].strip())
                    s4g_dict[galaxy]['flux_fracs'].append(flux_frac)

            if line[17:25].strip() == 'edgedisk':

                try:
                    if flux_frac > s4g_dict[galaxy]['flux_fracs'][-1]:
                        s4g_dict[galaxy]['h_r'] = float(float(line[52:59].strip()))
                        s4g_dict[galaxy]['flux_fracs'].append(flux_frac)
                except IndexError:
                    s4g_dict[galaxy]['h_r'] = float(float(line[52:59].strip()))
                    s4g_dict[galaxy]['flux_fracs'].append(flux_frac)

f.close()

final_galaxies = []
final_scale_length = []
final_bulge_total_ratio = []

# Now match these up with ALMA galaxies

for galaxy in alma_galaxies:

    galaxy = galaxy.strip()

    try:
        galaxy_dict = s4g_dict[galaxy]
    except KeyError:
        final_galaxies.append(galaxy)
        final_scale_length.append(np.nan)
        final_bulge_total_ratio.append(np.nan)
        continue

    final_galaxies.append(galaxy)

    try:
        final_scale_length.append(galaxy_dict['h_r'])
    except KeyError:
        final_scale_length.append(np.nan)

    try:
        final_bulge_total_ratio.append(galaxy_dict['b_t_ratio'])
    except KeyError:
        final_bulge_total_ratio.append(np.nan)

np.savetxt('lit_tables/s4g_scale_lengths.txt',
           np.c_[final_galaxies, final_scale_length, final_bulge_total_ratio],
           fmt='%s')

print('Complete!')
