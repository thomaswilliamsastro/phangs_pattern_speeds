# -*- coding: utf-8 -*-
"""
Test to see how varying how much effect the data reduction version has on the ALMA pattern speeds

@author: Tom Williams
"""

import os
import time

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from vars import phangs_folder, alma_plot, muse_plot, alma_output, muse_output, alma_galaxies, muse_galaxies, \
    hdu_types, star_masks, mask_outside_bars

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

start = time.time()

instrument = 'muse'

if instrument == 'alma':
    output_folder = alma_output
    plot_folder = alma_plot
    galaxies = alma_galaxies
    old_version = 'v33'
    new_version = 'v34'
elif instrument == 'muse':
    output_folder = muse_output
    plot_folder = muse_plot
    galaxies = muse_galaxies
    old_version = 'DR1.0'
    new_version = 'DR2.0'

hdu_type = hdu_types[0]
star_mask = star_masks[0]
mask_outside_bar = mask_outside_bars[0]

if not os.path.exists(plot_folder):
    os.mkdir(plot_folder)

old_mask_omegas = []
old_mask_omegas_up = []
old_mask_omegas_down = []

new_mask_omegas = []
new_mask_omegas_up = []
new_mask_omegas_down = []

for galaxy in galaxies:

    galaxy = galaxy.strip()

    file_name_old = output_folder + old_version + '/'

    if instrument == 'muse':
        file_name_old += hdu_type + '_'
        if star_mask:
            file_name_old += 'smask_'
        else:
            file_name_old += 'nosmask_'

    if mask_outside_bar:
        file_name_old += 'bmask/'
    else:
        file_name_old += 'nobmask/'

    file_name_old += galaxy + '_'

    if instrument == 'muse':
        file_name_old += hdu_type + '_'
        if star_mask:
            file_name_old += 'smask_'
        else:
            file_name_old += 'nosmask_'

    if mask_outside_bar:
        file_name_old += 'bmask_'
    else:
        file_name_old += 'nobmask_'

    file_name_old += 'pattern_speed_' + instrument + '.txt'
    file_name_new = file_name_old.replace(old_version, new_version)

    try:
        old_mask_omega, old_mask_omega_up, old_mask_omega_down = np.loadtxt(file_name_old, unpack=True)
        new_mask_omega, new_mask_omega_up, new_mask_omega_down = np.loadtxt(file_name_new, unpack=True)
        # print(galaxy, old_mask_omega/new_mask_omega)
    except OSError:
        print(galaxy)
        continue

    old_mask_omegas.append(old_mask_omega)
    old_mask_omegas_up.append(old_mask_omega_up)
    old_mask_omegas_down.append(old_mask_omega_down)

    new_mask_omegas.append(new_mask_omega)
    new_mask_omegas_up.append(new_mask_omega_up)
    new_mask_omegas_down.append(new_mask_omega_down)

plt.figure(figsize=(8, 6))
ax = plt.subplot(1, 1, 1)

# plt.subplot(2, 2, 1)

plt.errorbar(np.asarray(old_mask_omegas), np.asarray(new_mask_omegas),
             xerr=[np.asarray(old_mask_omegas_down), np.asarray(old_mask_omegas_up)],
             yerr=[np.asarray(new_mask_omegas_down), np.asarray(new_mask_omegas_up)],
             c='k', ls='none', marker='o')

plt.plot([0, 100], [0, 100],
         c='k', ls='--')

plt.xlim([0, 100])
plt.ylim([0, 100])

plt.xlabel(r'$\Omega_\mathrm{p, ' + old_version + '}$ (km s$^{-1}$ kpc$^{-1}$)')
plt.ylabel(r'$\Omega_\mathrm{p, ' + new_version + '}$ (km s$^{-1}$ kpc$^{-1}$)')

plt.tight_layout()

ax.set_aspect('equal')

# plt.show()

plt.savefig(plot_folder + 'reduction_comparison.png',
            bbox_inches='tight')
plt.savefig(plot_folder + 'reduction_comparison.pdf',
            bbox_inches='tight')

print('Complete! Took %.2fs' % (time.time() - start))
