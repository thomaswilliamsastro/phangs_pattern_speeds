# -*- coding: utf-8 -*-
"""
Compare pattern speeds versus corrected, with slits removed

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from astropy.table import Table

from vars import phangs_folder, alma_version, alma_output, alma_galaxies, alma_plot, plot_folder

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

use_lowres = False

galaxies = []

if use_lowres:
    extension = '_bmask_lowres'
else:
    extension = '_bmask'

pattern_speeds_dict = {}

null_names = ['moderate', 'stricter', 'strictest']
type_names = ['corr', 'null']

pattern_speeds_dict['original'] = np.array([])
pattern_speeds_dict['original_err_up'] = np.array([])
pattern_speeds_dict['original_err_down'] = np.array([])

for null_name in null_names:
    for type_name in type_names:
        pattern_speeds_dict[type_name + null_name] = np.array([])
        pattern_speeds_dict[type_name + null_name + '_err_up'] = np.array([])
        pattern_speeds_dict[type_name + null_name + '_err_down'] = np.array([])

for galaxy in alma_galaxies:

    galaxy = galaxy.upper()
    galaxies.append(galaxy)

    file_name = os.path.join(alma_output, alma_version, 'bmask', galaxy + extension + '_pattern_speed_alma.txt')
    if not os.path.exists(file_name):
        pattern_speed = [np.nan, np.nan, np.nan]
    else:
        pattern_speed = np.loadtxt(file_name)

    pattern_speeds_dict['original'] = np.append(pattern_speeds_dict['original'], pattern_speed[0])
    pattern_speeds_dict['original_err_up'] = np.append(pattern_speeds_dict['original_err_up'], pattern_speed[1])
    pattern_speeds_dict['original_err_down'] = np.append(pattern_speeds_dict['original_err_down'], pattern_speed[2])

    for null_name in null_names:

        for type_name in type_names:

            file_name = os.path.join(alma_output, alma_version, 'bmask', type_name, null_name,
                                     galaxy + extension + '_' + type_name)
            if type_name == 'corr':
                file_name += '_3_sigma'
            file_name += '_pattern_speed_alma.txt'

            if not os.path.exists(file_name):
                pattern_speed = [np.nan, np.nan, np.nan]
            else:
                pattern_speed = np.loadtxt(file_name)

            pattern_speeds_dict[type_name + null_name] = np.append(pattern_speeds_dict[type_name + null_name],
                                                                   pattern_speed[0])
            pattern_speeds_dict[type_name + null_name + '_err_up'] = np.append(
                pattern_speeds_dict[type_name + null_name + '_err_up'],
                pattern_speed[1])
            pattern_speeds_dict[type_name + null_name + '_err_down'] = np.append(
                pattern_speeds_dict[type_name + null_name + '_err_down'],
                pattern_speed[2])

# Save these out

table_file_name = os.path.join(alma_output, 'pattern_speeds_null_hypotheses.fits')

tab = Table()
tab.add_column(galaxies, name='galaxy')

for null_name in null_names:
    for type_name in type_names:

        tab.add_column(pattern_speeds_dict[type_name + null_name], name='_'.join([type_name, null_name]))
        tab.add_column(pattern_speeds_dict[type_name + null_name + '_err_up'],
                       name='_'.join([type_name, null_name, 'err_up']))
        tab.add_column(pattern_speeds_dict[type_name + null_name + '_err_down'],
                       name='_'.join([type_name, null_name, 'err_down']))

tab.write(table_file_name, overwrite=True)
no

pdf_file = plot_folder + 'diagnostics/null_comparisons'
if use_lowres:
    pdf_file += '_lowres'
pdf_file += '.pdf'

# with PdfPages(pdf_file) as pdf:
#     for galaxy in galaxies:
#
#         plt.figure(figsize=(12, 6))
#         plt.suptitle(galaxy)
#
#         for i, null_name in enumerate(null_names):
#
#             plot_name = os.path.join(alma_plot, alma_version, 'bmask', 'corr', null_name,
#                                      galaxy + extension + '_corr_3_sigma_alma.png')
#
#             if os.path.exists(plot_name):
#                 plt.subplot(1, 3, i + 1)
#                 plt.title(null_name)
#                 plt.imshow(plt.imread(plot_name))
#                 plt.axis('off')
#
#         plt.tight_layout()
#         # plt.show()
#         pdf.savefig(dpi=300)
#         plt.close()

# galaxies = np.array(galaxies)
# print(median_absolute_deviation((pattern_speeds - pattern_speeds_corr) / pattern_speeds))
#
# different_idx = np.where(np.abs(pattern_speeds - pattern_speeds_corr) >=
#                          np.sqrt(pattern_speeds_err_down ** 2 + pattern_speeds_corr_err_down ** 2))
# print(galaxies[different_idx])

xlims = [0, 100]
c = ['k', 'g', 'r']

plot_name = os.path.join(alma_plot, alma_version, 'alma_comparison_corr')
if use_lowres:
    plot_name += '_lowres'

plt.figure(figsize=(6, 4))

for i, null_name in enumerate(null_names):

    plt.errorbar(
        pattern_speeds_dict['original'], pattern_speeds_dict['corr' + null_name],
        xerr=[
            pattern_speeds_dict['original_err_down'], pattern_speeds_dict['original_err_up']
        ],
        yerr=[
            pattern_speeds_dict['corr' + null_name + '_err_down'], pattern_speeds_dict['corr' + null_name + '_err_up']
        ],
        c=c[i], ls='none', marker='o', label=null_name)

plt.plot(xlims, xlims, c='k', ls='--')

plt.xlim(xlims)
plt.ylim(xlims)

plt.xlabel(r'$\Omega_\mathrm{p, corr}$ (km s$^{-1}$ kpc$^{-1}$)')
plt.ylabel(r'$\Omega_\mathrm{p, null}$ (km s$^{-1}$ kpc$^{-1}$)')

plt.legend(loc='upper left', frameon=False)

plt.tight_layout()
# plt.show()

plt.savefig(plot_name + '.png', bbox_inches='tight')
plt.savefig(plot_name + '.pdf', bbox_inches='tight')

plt.close()

# And do a lil KDE plot of the difference

plot_name = os.path.join(alma_plot, alma_version, 'alma_comparison_corr_kde')
if use_lowres:
    plot_name += '_lowres'

plt.figure(figsize=(4, 3))

for i, null_name in enumerate(null_names):

    pattern_speeds_residual = \
        (pattern_speeds_dict['corr' + null_name] - pattern_speeds_dict['original']) / pattern_speeds_dict['original']
    pattern_speeds_residual[np.abs(pattern_speeds_residual) > 2] = np.nan

    sns.kdeplot(pattern_speeds_residual, color=c[i], bw='silverman', label=null_name)

plt.legend(loc='upper left')

plt.xlabel(r'$(\Omega_\mathrm{p, corr} - \Omega_p) / \Omega_p$')
plt.ylabel('Probability Density')

# plt.show()

plt.tight_layout()

plt.savefig(plot_name + '.png', bbox_inches='tight')
plt.savefig(plot_name + '.pdf', bbox_inches='tight')

plt.close()

# Finally, look at the difference between the corrected values and the values from the nulls, for each null.

plot_name = os.path.join(alma_plot, alma_version, 'alma_comparison_corr_to_null')
if use_lowres:
    plot_name += '_lowres'

plt.figure(figsize=(6, 4))

for i, null_name in enumerate(null_names):

    plt.errorbar(
        pattern_speeds_dict['corr' + null_name], pattern_speeds_dict['null' + null_name],
        xerr=[
            pattern_speeds_dict['corr' + null_name + '_err_down'], pattern_speeds_dict['corr' + null_name + '_err_up']
        ],
        yerr=[
            pattern_speeds_dict['null' + null_name + '_err_down'], pattern_speeds_dict['null' + null_name + '_err_up']
        ],
        c=c[i], ls='none', marker='o', label=null_name)

plt.plot(xlims, xlims, c='k', ls='--')

plt.xlim(xlims)
plt.ylim(xlims)

plt.xlabel(r'$\Omega_\mathrm{p, corr}$ (km s$^{-1}$ kpc$^{-1}$)')
plt.ylabel(r'$\Omega_\mathrm{p, null}$ (km s$^{-1}$ kpc$^{-1}$)')

plt.legend(loc='upper left', frameon=False)

plt.tight_layout()

# plt.show()

plt.savefig(plot_name + '.png', bbox_inches='tight')
plt.savefig(plot_name + '.pdf', bbox_inches='tight')

plt.close()

print('Complete!')
