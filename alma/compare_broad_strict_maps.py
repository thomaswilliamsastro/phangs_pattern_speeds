# -*- coding: utf-8 -*-
"""
Compare broad and strict masks

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from astropy.table import Table

from vars import phangs_folder, alma_version, alma_plot, alma_output, alma_galaxies, mask_outside_bars, output_folder, \
    corot_version

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

if not os.path.exists(alma_plot):
    os.mkdir(alma_plot)
if not os.path.exists(alma_plot + alma_version):
    os.mkdir(alma_plot + alma_version)

mask_outside_bar = mask_outside_bars[0]

omega_bars_residuals = []
galaxies = []
q_flags = []

corot_t = Table.read(output_folder + 'pattern_speed_table_' + corot_version + '.fits')

for galaxy in alma_galaxies:

    file_name = alma_output + alma_version

    if mask_outside_bar:
        file_name += '/bmask/'
    else:
        file_name += '/nobmask/'

    galaxy = galaxy.strip()

    q_flag = corot_t[corot_t['GALAXY'] == galaxy]['OM_P_ALMA_QUAL']

    file_name += galaxy

    if mask_outside_bar:
        file_name += '_bmask'
    else:
        file_name += '_nobmask'

    file_name += '_pattern_speed_alma.txt'

    try:
        omega_bar_broad = np.loadtxt(file_name, usecols=0)
        omega_bar_strict = np.loadtxt(file_name.replace('mask_', 'mask_strict_'), usecols=0)

        omega_bars_residuals.append((omega_bar_broad - omega_bar_strict) / omega_bar_broad)
        galaxies.append(galaxy)
        try:
            q_flags.append(q_flag[0])
        except IndexError:
            q_flags.append(3)
    except OSError:
        continue

# Mask any exceptionally deviant differences

for i, galaxy in enumerate(galaxies):

    print(galaxy, omega_bars_residuals[i])

no

omega_bars_residuals = np.array(omega_bars_residuals)
q_flags = np.array(q_flags)

mask = np.where(np.abs(omega_bars_residuals) < 1)
omega_bars_residuals = omega_bars_residuals[mask]
q_flags = q_flags[mask]

q_flag_mask = np.where(q_flags <= 2)

plt.figure(figsize=(8, 6))

sns.kdeplot(omega_bars_residuals, bw='silverman', color='k', shade=False, label='All')
sns.kdeplot(omega_bars_residuals[q_flag_mask], bw='silverman', color='b', shade=True, label='$Q = 1,2$')

plt.ylabel('Probability Density')
plt.xlabel(r'$\frac{\Omega_\mathrm{p, broad} - \Omega_\mathrm{p, strict}}{\Omega_\mathrm{p, broad}}$')

plt.legend(loc='upper right')

plot_name = alma_plot + alma_version + '/broad_strict_comparison'

plt.savefig(plot_name + '.png', bbox_inches='tight')
plt.savefig(plot_name + '.pdf', bbox_inches='tight')
plt.close()

print('Complete!')
