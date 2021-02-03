# -*- coding: utf-8 -*-
"""
Plot distributions for the streaming motion mocks

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.pyplot import cm

from vars import phangs_folder, output_folder, plot_folder

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

use_full_cov = True
compare_to_true = True

plot_filename = plot_folder + 'oct2020mocks'
if use_full_cov:
    plot_filename += '_fullcov'
    title = 'Fully Covered'
else:
    title = 'Incomplete Coverage'

plot_filename += '/streaming_motions_comparisons'
if compare_to_true:
    plot_filename += '_true_comp'

bar_lens = ['10.0', '30.0', '50.0']
bar_ellips = ['0.25', '0.5']
bar_angs = ['60.0', '90.0', '120.0', '150.0', '180.0', '210.0']
fs = ['0.0', '0.25', '0.5', '0.75', '1.0']

n = len(bar_lens) * len(fs) * len(bar_ellips)

colours = iter(cm.rainbow(np.linspace(0, 1, n)))

plt.figure(figsize=(4, 3))

for f in fs:
    for bar_len in bar_lens:
        for bar_ellip in bar_ellips:

            om_clumps = []

            for bar_ang in bar_angs:
                file_name = output_folder + 'oct2020mocks'
                if use_full_cov:
                    file_name += '_fullcov'

                file_name += '/' + bar_len + '_' + bar_ellip + '_' + bar_ang + '_' + f + '_pattern_speed.txt'
                if not os.path.exists(file_name):
                    continue

                om_clump = np.loadtxt(file_name, usecols=0)
                om_clumps.append(om_clump)

            om_clumps = np.array(om_clumps)

            if compare_to_true:

                if bar_len == '10.0':
                    om_clumps_percent = om_clumps / (76.2 * 9.84/10)
                if bar_len == '30.0':
                    om_clumps_percent = om_clumps / (61.5 * 9.84/10)
                if bar_len == '50.0':
                    om_clumps_percent = om_clumps / (48.2 * 9.84/10)

            else:

                om_clumps_percent = om_clumps / np.nanmean(om_clumps)

            colour = next(colours)
            sns.kdeplot(om_clumps_percent, color=colour, shade=False, bw='silverman',
                        #label='%s, %s' % (f, bar_len),
                        )
            print(f, bar_len, np.nanpercentile(om_clumps_percent, 84) - np.nanmedian(om_clumps_percent))

# plt.legend(loc='upper right')

plt.xlim([0.2, 0.8])

plt.xlabel(r'$\Omega_p/\overline{\Omega_p}$')
plt.ylabel('Probability Density')

plt.title(title)

plt.tight_layout()

# plt.show()

plt.savefig(plot_filename + '.png', bbox_inches='tight')
plt.savefig(plot_filename + '.pdf', bbox_inches='tight')

plt.close()

print('Complete!')
