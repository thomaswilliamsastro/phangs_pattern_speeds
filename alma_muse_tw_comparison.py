# -*- coding: utf-8 -*-
"""
Compare MUSE/ALMA TW pattern speeds

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from astropy.table import Table
from scipy.stats import median_absolute_deviation

from vars import phangs_folder, alma_version, muse_version, muse_galaxies, plot_folder, output_folder, corot_version

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

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

galaxies_mad = []

good_alma_rows = corot_table[(corot_table['OM_P_ALMA_QUAL'] == 1) | (corot_table['OM_P_ALMA_QUAL'] == 2)]

alma_om_ps = []
alma_qs = []
muse_ha_om_ps = []
muse_ha_qs = []
muse_mass_om_ps = []
muse_mass_qs = []

for galaxy in corot_table['GALAXY']:

    galaxy_mad = []

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
    alma_om_ps.append(alma_om)
    muse_ha_om_ps.append(muse_om_ha)
    muse_mass_om_ps.append(muse_om_mass)

    alma_qs.append(alma_om_q)
    muse_ha_qs.append(muse_om_ha_q)
    muse_mass_qs.append(muse_om_mass_q)

    if muse_om_mass_q in [1, 2]:
        fill_style = 'full'
        galaxy_mad.append(muse_om_mass)
    else:
        fill_style = 'none'

    plt.errorbar(muse_om_mass, position - 0.125,
                 xerr=np.array([[muse_om_mass_err_down, muse_om_mass_err_up]]).T,
                 fmt='o', c='r', fillstyle=fill_style)

    if muse_om_ha_q in [1, 2]:
        fill_style = 'full'
        galaxy_mad.append(muse_om_ha)
    else:
        fill_style = 'none'

    plt.errorbar(muse_om_ha, position,
                 xerr=np.array([[muse_om_ha_err_down, muse_om_ha_err_up]]).T,
                 fmt='o', c='cyan', fillstyle=fill_style)

    if alma_om_q in [1, 2]:
        fill_style = 'full'
        galaxy_mad.append(alma_om)
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

    if len(galaxy_mad) > 0:
        galaxy_mad = np.array(galaxy_mad)
        galaxy_mean = np.mean(galaxy_mad)
        galaxy_mad -= galaxy_mean
        galaxy_mad /= galaxy_mean

        galaxies_mad.extend(galaxy_mad)

alma_om_ps = np.array(alma_om_ps)
muse_ha_om_ps = np.array(muse_ha_om_ps)
muse_mass_om_ps = np.array(muse_mass_om_ps)
alma_qs = np.array(alma_qs)
muse_ha_qs = np.array(muse_ha_qs)
muse_mass_qs = np.array(muse_mass_qs)
galaxies_mad = np.array(galaxies_mad)

print('Percentage MAD: %.2f' % (median_absolute_deviation(galaxies_mad) * 100))

plt.yticks(0.5 * np.array(range(len(galaxies))), galaxies)
plt.xlabel(r'$\Omega_{p}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')

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

# Look at the residual differences between these things

plt.figure(figsize=(8, 4))

plt.subplot(1, 3, 1)

alma_muse_ha_residuals = (alma_om_ps - muse_ha_om_ps) / alma_om_ps
q_flags = ((alma_qs < 3) & (muse_ha_qs < 3))
mask = np.where(np.abs(alma_muse_ha_residuals) < 1)
omega_bars_residuals = alma_muse_ha_residuals[mask]
q_flags = q_flags[mask]

sns.kdeplot(omega_bars_residuals, color='k', shade=False, bw='silverman')
sns.kdeplot(omega_bars_residuals[q_flags], color='r', shade=True, bw='silverman')

plt.xlim([-1, 1])

plt.ylabel('Probability Density')
plt.xlabel(r'$\frac{\Omega_\mathrm{p, ALMA} - \Omega_\mathrm{p, MUSE\,H\alpha}}{\Omega_\mathrm{p, ALMA}}$')

plt.subplot(1, 3, 2)

alma_muse_mass_residuals = (alma_om_ps - muse_mass_om_ps) / alma_om_ps
q_flags = ((alma_qs < 3) & (muse_mass_qs < 3))
mask = np.where(np.abs(alma_muse_mass_residuals) < 1)
omega_bars_residuals = alma_muse_mass_residuals[mask]
q_flags = q_flags[mask]

sns.kdeplot(omega_bars_residuals, color='k', shade=False, bw='silverman')
sns.kdeplot(omega_bars_residuals[q_flags], color='r', shade=True, bw='silverman')

plt.xlim([-1, 1])

plt.xlabel(r'$\frac{\Omega_\mathrm{p, ALMA} - \Omega_\mathrm{p, MUSE\,M_\star}}{\Omega_\mathrm{p, ALMA}}$')

plt.subplot(1, 3, 3)

muse_ha_muse_mass_residuals = (muse_ha_om_ps - muse_mass_om_ps) / muse_ha_om_ps
q_flags = ((muse_ha_qs < 3) & (muse_mass_qs < 3))
mask = np.where(np.abs(muse_ha_muse_mass_residuals) < 1)
omega_bars_residuals = muse_ha_muse_mass_residuals[mask]
q_flags = q_flags[mask]

sns.kdeplot(omega_bars_residuals, color='k', shade=False, bw='silverman', label='All')
sns.kdeplot(omega_bars_residuals[q_flags], color='r', shade=True, bw='silverman', label=r'$Q=1,2$')

plt.xlim([-1, 1])

plt.xlabel(r'$\frac{\Omega_\mathrm{p, MUSE\,H\alpha} - \Omega_\mathrm{p, MUSE\,M_\star}}'
           r'{\Omega_\mathrm{p, MUSE\,H\alpha}}$')

plt.legend(loc='upper right')

plt.tight_layout()

plot_filename += '_kde'

plt.savefig(plot_filename + '.png',
            bbox_inches='tight')
plt.savefig(plot_filename + '.pdf',
            bbox_inches='tight')

print('Complete!')
