# -*- coding: utf-8 -*-
"""
Check effects of FOV on pattern speed

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from vars import phangs_folder, muse_version, muse_output, muse_galaxies, plot_folder

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

galaxies = []

pattern_speeds = np.array([])
pattern_speeds_err_up = np.array([])
pattern_speeds_err_down = np.array([])

for galaxy in muse_galaxies:

    pattern_speed_filename = os.path.join(muse_output, muse_version, 'mass_smask_bmask', 'fov_check',
                                          galaxy + '_mass_smask_bmask_fov_check_pattern_speed_muse.txt')

    if not os.path.exists(pattern_speed_filename):
        continue

    pattern_speed = np.loadtxt(pattern_speed_filename)

    galaxies.append(galaxy.upper())

    pattern_speeds = np.append(pattern_speeds, pattern_speed[0])
    pattern_speeds_err_up = np.append(pattern_speeds_err_up, pattern_speed[1])
    pattern_speeds_err_down = np.append(pattern_speeds_err_down, pattern_speed[2])

x = range(len(galaxies))

plot_name = os.path.join(plot_folder, 'muse', muse_version, 'muse_fov_comparison')

plt.figure(figsize=(6, 4))

plt.errorbar(x, pattern_speeds, yerr=[pattern_speeds_err_down, pattern_speeds_err_up],
             c='k', ls='none', marker='o')

plt.xticks(x, galaxies, rotation=45)
plt.ylabel(r'$\Omega_\mathrm{p}$ (km s$^{-1}$ kpc$^{-1}$)')

plt.tight_layout()

# plt.show()

plt.savefig(plot_name + '.png', bbox_inches='tight')
plt.savefig(plot_name + '.pdf', bbox_inches='tight')

plt.close()

print('Complete!')
