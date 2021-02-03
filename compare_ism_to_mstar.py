# -*- coding: utf-8 -*-
"""
Published literature pattern speeds comparison to this work

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from astropy.table import Table

from vars import phangs_folder, output_folder, pattern_speed_version, plot_folder

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

pattern_speeds_table = Table.read(output_folder + 'pattern_speed_table_' + pattern_speed_version + '.fits')

co_rows = pattern_speeds_table[(
        ((pattern_speeds_table['OM_P_MUSE_MASS_QUAL'] == 1) | (pattern_speeds_table['OM_P_MUSE_MASS_QUAL'] == 2)) &
        ((pattern_speeds_table['OM_P_ALMA_QUAL'] == 1) | (pattern_speeds_table['OM_P_ALMA_QUAL'] == 2))
)]
co_diff = (co_rows['OM_P_ALMA'] - co_rows['OM_P_MUSE_MASS']) / co_rows['OM_P_MUSE_MASS']

ha_rows = pattern_speeds_table[(
        ((pattern_speeds_table['OM_P_MUSE_HA_QUAL'] == 1) | (pattern_speeds_table['OM_P_MUSE_HA_QUAL'] == 2)) &
        ((pattern_speeds_table['OM_P_MUSE_MASS_QUAL'] == 1) | (pattern_speeds_table['OM_P_MUSE_MASS_QUAL'] == 2))
)]
ha_diff = (ha_rows['OM_P_MUSE_HA'] - ha_rows['OM_P_MUSE_MASS']) / ha_rows['OM_P_MUSE_MASS']

median_co_diff = np.nanmedian(co_diff)
median_ha_diff = np.nanmedian(ha_diff)

print(median_co_diff, median_ha_diff)

plot_name = plot_folder + 'ism_mstar_om_p_comparison'

plt.figure(figsize=(4, 3))

sns.kdeplot(co_diff, color='b', shade=False, label='CO, N=%d' % len(co_rows))
sns.kdeplot(ha_diff, color='cyan', shade=False, label=r'H$\alpha$, N=%d' % len(ha_rows))

plt.axvline(0, c='k', ls='--')
plt.axvline(median_co_diff, c='b', ls='--')
plt.axvline(median_ha_diff, c='cyan', ls='--')

plt.xlabel(r'Pattern Speed Difference')
plt.ylabel('Probability Density')

plt.xticks([-0.5, 0, 0.5, 1, 1.5])

plt.legend(frameon=False)

plt.tight_layout()
# plt.show()

plt.savefig(plot_name + '.png', bbox_inches='tight')
plt.savefig(plot_name + '.pdf', bbox_inches='tight')

plt.close()

print('Complete!')
