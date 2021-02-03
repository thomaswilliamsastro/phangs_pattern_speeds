# -*- coding: utf-8 -*-
"""
Plot the ALMA pattern speeds for all the PHANGS galaxies

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from astropy.table import Table

from vars import phangs_folder, alma_version, alma_plot, pattern_speed_version, output_folder, plot_folder

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

if not os.path.exists(alma_plot):
    os.mkdir(alma_plot)
if not os.path.exists(alma_plot + alma_version):
    os.mkdir(alma_plot + alma_version)

# Plot including quality flag information from the pattern speed table

use_muse = True

corot_table = Table.read(output_folder + 'pattern_speed_table_' + pattern_speed_version + '.fits')

if use_muse:
    good_alma_om_p = corot_table[(corot_table['OM_P_MUSE_MASS_QUAL'] == 1) |
                                 (corot_table['OM_P_MUSE_MASS_QUAL'] == 2)]['OM_P_MUSE_MASS']
    all_alma_om_p = corot_table['OM_P_MUSE_MASS']
    plot_name = plot_folder + 'muse_om_p_distribution'
else:
    good_alma_om_p = corot_table[(corot_table['OM_P_ALMA_QUAL'] == 1) |
                                 (corot_table['OM_P_ALMA_QUAL'] == 2)]['OM_P_ALMA']
    all_alma_om_p = corot_table['OM_P_ALMA']
    plot_name = plot_folder + 'alma_om_p_distribution'

print(np.nanpercentile(good_alma_om_p, [16, 50, 84]))
print(np.nanpercentile(all_alma_om_p, [16, 50, 84]))

plt.figure(figsize=(4, 3))
sns.kdeplot(all_alma_om_p,
            bw='silverman', color='r', shade=False,
            clip=[0, 100], label='All')
sns.kdeplot(good_alma_om_p,
            bw='silverman', color='k', shade=True,
            clip=[0, 100], label='$Q = 1,2$')

plt.ylabel('Probability Density')
plt.xlabel(r'$\Omega_\mathrm{P}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')

plt.savefig(plot_name + '.png',
            bbox_inches='tight')
plt.savefig(plot_name + '.pdf',
            bbox_inches='tight')
plt.close()

print('Complete!')
