# -*- coding: utf-8 -*-
"""
Compare the values calculated from using MUSE stellar mass vs. IRAC 3.6micron image

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from vars import phangs_folder, muse_version, muse_plot, muse_output, muse_galaxies

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

if not os.path.exists(muse_plot):
    os.mkdir(muse_plot)
if not os.path.exists(muse_plot + muse_version):
    os.mkdir(muse_plot + muse_version)

mask_outside_bar = True
star_mask = True
hdu_type = 'mass'

muse_omegas = []
muse_omegas_up = []
muse_omegas_down = []

irac_omegas = []
irac_omegas_up = []
irac_omegas_down = []

for galaxy in muse_galaxies:

    galaxy = galaxy.strip()

    file_name_muse = muse_output + muse_version + '/'

    file_name_muse += hdu_type + '_'
    if star_mask:
        file_name_muse += 'smask_'
    else:
        file_name_muse += 'nosmask_'

    if mask_outside_bar:
        file_name_muse += 'bmask/'
    else:
        file_name_muse += 'nobmask/'

    file_name_muse += galaxy + '_'

    file_name_muse += hdu_type + '_'
    if star_mask:
        file_name_muse += 'smask_'
    else:
        file_name_muse += 'nosmask_'

    if mask_outside_bar:
        file_name_muse += 'bmask_'
    else:
        file_name_muse += 'nobmask_'

    file_name_muse += 'pattern_speed_muse.txt'
    file_name_irac = file_name_muse.replace('mass', 'spitzer_mass')

    try:
        muse_omega, muse_omega_up, muse_omega_down = np.loadtxt(file_name_muse, unpack=True)
        irac_omega, irac_omega_up, irac_omega_down = np.loadtxt(file_name_irac, unpack=True)
        # print(galaxy, old_mask_omega/new_mask_omega)
    except OSError:
        print(galaxy)
        continue

    muse_omegas.append(muse_omega)
    muse_omegas_up.append(muse_omega_up)
    muse_omegas_down.append(muse_omega_down)

    irac_omegas.append(irac_omega)
    irac_omegas_up.append(irac_omega_up)
    irac_omegas_down.append(irac_omega_down)

plt.figure(figsize=(8, 6))
ax = plt.subplot(1, 1, 1)

# plt.subplot(2, 2, 1)

plt.errorbar(np.asarray(muse_omegas), np.asarray(irac_omegas),
             xerr=[np.asarray(muse_omegas_down), np.asarray(muse_omegas_up)],
             yerr=[np.asarray(irac_omegas_down), np.asarray(irac_omegas_up)],
             c='k', ls='none', marker='o')

plt.plot([0, 100], [0, 100],
         c='k', ls='--')

plt.xlim([0, 100])
plt.ylim([0, 100])

plt.xlabel(r'$\Omega_\mathrm{p, MUSE}$ (km s$^{-1}$ kpc$^{-1}$)')
plt.ylabel(r'$\Omega_\mathrm{p, IRAC}$ (km s$^{-1}$ kpc$^{-1}$)')

plt.tight_layout()

ax.set_aspect('equal')

plt.show()

print('Complete!')
