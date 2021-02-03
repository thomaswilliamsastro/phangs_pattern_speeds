# -*- coding: utf-8 -*-
"""
Compare the various MUSE emission options to see if there's significant differences

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
if not os.path.exists(muse_plot + muse_version + '/emission/'):
    os.mkdir(muse_plot + muse_version + '/emission/')

labels = [r'$M_\ast$', 'Flux', r'H$\alpha$']

methods = []

for hdu_type in ['mass', 'flux', 'ha']:

    for smask in [False, True]:

        for bmask in [False, True]:

            if smask:
                star_ext = '_smask'
            else:
                star_ext = '_nosmask'

            if bmask:
                mask_ext = '_bmask'
            else:
                mask_ext = '_nobmask'

            methods.append(hdu_type + star_ext + mask_ext)

for galaxy in muse_galaxies:

    galaxy_found = True

    print('Starting %s' % galaxy)

    plt.figure(figsize=(8, 6))
    ax = plt.subplot(111)

    position = 0

    for method in methods:

        file_name = muse_output + muse_version + '/' + method + '/' + galaxy + '_' + method + '_pattern_speed_muse.txt'

        try:
            omega_bar, omega_bar_err_up, omega_bar_err_down = np.loadtxt(file_name, unpack=True)
        except OSError:
            galaxy_found = False
            continue

        if '_nosmask_nobmask' in method:
            c = 'r'

            plt.axhline(position + 0.375,
                        c='k',
                        ls='--')

            if method == methods[0]:
                plt.errorbar(-100, -100,
                             xerr=1,
                             fmt='o',
                             c=c,
                             label='No star mask, no bar mask')

        elif '_nosmask_bmask' in method:
            c = 'g'

            if method == methods[1]:
                plt.errorbar(-100, -100,
                             xerr=1,
                             fmt='o',
                             c=c,
                             label='No star mask, bar mask')

        elif '_smask_nobmask' in method:
            c = 'gold'

            if method == methods[2]:
                plt.errorbar(-100, -100,
                             xerr=1,
                             fmt='o',
                             c=c,
                             label='Star mask, no bar mask')

        else:
            c = 'b'

            if method == methods[3]:
                plt.errorbar(-100, -100,
                             xerr=1,
                             fmt='o',
                             c=c,
                             label='Star mask, bar mask')

        plt.errorbar(omega_bar, position - 0.125 / 2,
                     xerr=np.array([[omega_bar_err_down, omega_bar_err_up]]).T,
                     fmt='o',
                     c=c)

        position += 0.125

    if not galaxy_found:
        plt.close()
        print('%s not found: skipping' % galaxy)
        continue

    plt.ylim([-0.125, 0.5 * len(labels) + 0.5 - 0.125 / 2])

    xlims = plt.xlim()
    plt.xlim([0.4 * np.nanmin(omega_bar), xlims[-1]])

    plt.yticks(0.5 * np.array(range(len(labels))) + 0.125, labels)
    plt.xlabel(r'$\Omega_{p}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')

    plt.text(0.05, 0.95, galaxy,
             ha='left', va='top',
             fontweight='bold',
             transform=ax.transAxes)

    plt.legend(loc='upper right',
               frameon=False)

    plt.tight_layout()

    plot_file_name = muse_plot + muse_version + '/emission/' + galaxy + '_muse_emission_comparison'

    plt.savefig(plot_file_name + '.png',
                bbox_inches='tight')
    plt.savefig(plot_file_name + '.pdf',
                bbox_inches='tight')

    plt.close()

print('Complete!')
