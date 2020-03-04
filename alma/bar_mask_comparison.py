# -*- coding: utf-8 -*-
"""
Compare the effects of bar masking on the measured pattern speeds

@author: Tom Williams
"""

import os
import time

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from vars import phangs_folder, alma_version, alma_output, alma_plot, alma_galaxies

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

start = time.time()

os.chdir(phangs_folder)

if not os.path.exists(alma_plot):
    os.mkdir(alma_plot)
if not os.path.exists(alma_plot + alma_version):
    os.mkdir(alma_plot + alma_version)

mask_omegas = []
mask_omegas_up = []
mask_omegas_down = []

unmask_omegas = []
unmask_omegas_up = []
unmask_omegas_down = []

for galaxy in alma_galaxies:
    galaxy = galaxy.strip()

    try:
        mask_omega, mask_omega_up, mask_omega_down = np.loadtxt(
            alma_output + alma_version + '/bmask/' + galaxy + '_bmask_pattern_speed_alma.txt',
            unpack=True)
        unmask_omega, unmask_omega_up, unmask_omega_down = np.loadtxt(
            alma_output + alma_version + '/nobmask/' + galaxy + '_nobmask_pattern_speed_alma.txt',
            unpack=True)
    except OSError:
        continue

    mask_omegas.append(mask_omega)
    mask_omegas_up.append(mask_omega_up)
    mask_omegas_down.append(mask_omega_down)

    unmask_omegas.append(unmask_omega)
    unmask_omegas_up.append(unmask_omega_up)
    unmask_omegas_down.append(unmask_omega_down)

plt.figure(figsize=(8, 6))

# plt.subplot(2, 2, 1)

plt.errorbar(np.asarray(mask_omegas), np.asarray(unmask_omegas),
             xerr=[np.asarray(mask_omegas_down), np.asarray(mask_omegas_up)],
             yerr=[np.asarray(unmask_omegas_down), np.asarray(unmask_omegas_up)],
             c='k', ls='none', marker='o')

plt.plot([0, 100], [0, 100],
         c='k', ls='--')

plt.xlim([0, 100])
plt.ylim([0, 100])

plt.xlabel(r'$\Omega_\mathrm{p, mask}$ (km s$^{-1}$ kpc$^{-1}$)')
plt.ylabel(r'$\Omega_\mathrm{p, unmask}$ (km s$^{-1}$ kpc$^{-1}$)')

plt.tight_layout()

plt.savefig(alma_plot + alma_version + '/bar_mask_comparison.png',
            bbox_inches='tight')
plt.savefig(alma_plot + alma_version + '/bar_mask_comparison.pdf',
            bbox_inches='tight')

print('Complete! Took %.2fs' % (time.time() - start))
