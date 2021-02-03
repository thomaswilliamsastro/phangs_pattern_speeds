# -*- coding: utf-8 -*-
"""
Compare pattern speeds from the actual data versus the axisymmetric rotation velocities

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

pattern_speeds_obs = np.array([])
pattern_speeds_obs_err_up = np.array([])
pattern_speeds_obs_err_down = np.array([])

pattern_speeds_model = np.array([])
pattern_speeds_model_err_up = np.array([])
pattern_speeds_model_err_down = np.array([])

for galaxy in muse_galaxies:

    galaxy = galaxy.upper()
    if galaxy == 'NGC0628':
        galaxy = 'NGC628'

    pattern_speed_filename = os.path.join(muse_output, muse_version, 'ha_smask_bmask',
                                          galaxy + '_ha_smask_bmask_pattern_speed_muse.txt')
    pattern_speed_filename_model = os.path.join(muse_output, muse_version, 'ha_smask_bmask', 'null',
                                                galaxy + '_ha_smask_bmask_null_pattern_speed_muse.txt')

    if not os.path.exists(pattern_speed_filename) or not os.path.exists(pattern_speed_filename_model):
        continue

    pattern_speed = np.loadtxt(pattern_speed_filename)
    pattern_speed_model = np.loadtxt(pattern_speed_filename_model)

    galaxies.append(galaxy)

    pattern_speeds_obs = np.append(pattern_speeds_obs, pattern_speed[0])
    pattern_speeds_obs_err_up = np.append(pattern_speeds_obs_err_up, pattern_speed[1])
    pattern_speeds_obs_err_down = np.append(pattern_speeds_obs_err_down, pattern_speed[2])

    pattern_speeds_model = np.append(pattern_speeds_model, pattern_speed_model[0])
    pattern_speeds_model_err_up = np.append(pattern_speeds_model_err_up, pattern_speed_model[1])
    pattern_speeds_model_err_down = np.append(pattern_speeds_model_err_down, pattern_speed_model[2])

# with PdfPages(plot_folder + 'diagnostics/comparison_to_null.pdf') as pdf:
#     for galaxy in galaxies:
#         real_plot_name = os.path.join(alma_plot, alma_version, 'bmask', galaxy + '_bmask_alma.png')
#         model_plot_name = os.path.join(alma_plot, alma_version, 'bmask', 'null', galaxy + '_bmask_null_alma.png')
#
#         plt.figure(figsize=(12, 6))
#         plt.suptitle(galaxy)
#
#         plt.subplot(1, 2, 1)
#
#         plt.title('Obs')
#         plt.imshow(plt.imread(real_plot_name))
#         plt.axis('off')
#
#         plt.subplot(1, 2, 2)
#
#         plt.title('Model')
#         plt.imshow(plt.imread(model_plot_name))
#         plt.axis('off')
#
#         plt.tight_layout()
#         # plt.show()
#         pdf.savefig(dpi=300)
#         plt.close()

galaxies = np.array(galaxies)
# print(galaxies)

significant_diff_idx = np.where(np.abs(pattern_speeds_obs - pattern_speeds_model) <=
                                np.sqrt(pattern_speeds_model_err_down ** 2 + pattern_speeds_obs_err_down ** 2))
print(', '.join(list(galaxies[significant_diff_idx])))

print(len(significant_diff_idx[0]) / len(galaxies))

xlims = [0, 100]

plot_name = os.path.join(plot_folder, 'muse', muse_version, 'muse_comparison_null_hypothesis')

plt.figure(figsize=(6, 4))
plt.errorbar(pattern_speeds_model, pattern_speeds_obs,
             xerr=[pattern_speeds_model_err_down, pattern_speeds_model_err_up],
             yerr=[pattern_speeds_obs_err_down, pattern_speeds_obs_err_up],
             c='k', ls='none', marker='o')
plt.plot(xlims, xlims, c='k', ls='--')

plt.xlim(xlims)
plt.ylim(xlims)

plt.xlabel(r'$\Omega_\mathrm{p, obs}$ (km s$^{-1}$ kpc$^{-1}$)')
plt.ylabel(r'$\Omega_\mathrm{p, null}$ (km s$^{-1}$ kpc$^{-1}$)')

plt.tight_layout()
plt.show()

plt.savefig(plot_name + '.png', bbox_inches='tight')
plt.savefig(plot_name + '.pdf', bbox_inches='tight')

plt.close()

print('Complete!')
