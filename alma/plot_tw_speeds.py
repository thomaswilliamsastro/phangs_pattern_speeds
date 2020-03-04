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
from matplotlib import cm

from vars import phangs_folder, alma_version, alma_plot, alma_output, alma_galaxies, use_sim, mask_outside_bars

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

if not os.path.exists(alma_plot):
    os.mkdir(alma_plot)
if not os.path.exists(alma_plot + alma_version):
    os.mkdir(alma_plot + alma_version)

for mask_outside_bar in mask_outside_bars:

    output_folder = alma_output + alma_version + '/'

    if mask_outside_bar:
        output_folder += 'bmask/'
    else:
        output_folder += 'nobmask/'

    if use_sim:
        output_folder += 'sim/'

    plt.figure(figsize=(4, 0.5 * len(alma_galaxies)))
    ax1 = plt.subplot(1, 1, 1)

    frame1 = plt.gca()

    position = 0

    colours = iter(cm.rainbow(np.linspace(0, 1, len(alma_galaxies))))

    omega_bars = []

    for galaxy in alma_galaxies:

        galaxy = galaxy.strip()

        c = next(colours)

        if mask_outside_bar:
            file_name = galaxy + '_bmask'
        else:
            file_name = galaxy + '_nobmask'

        file_name += '_pattern_speed_alma.txt'

        try:
            omega_bar, omega_bar_err_up, omega_bar_err_down = np.loadtxt(output_folder + file_name,
                                                                         unpack=True)

            if use_sim:
                omega_bar -= 50

            omega_bars.append(omega_bar)

            plt.errorbar(omega_bar, position,
                         xerr=np.array([[omega_bar_err_down, omega_bar_err_up]]).T,
                         fmt='o',
                         c=c)
        except OSError:

            print(galaxy + ' not found')
            pass

        position += 0.25

    if use_sim:
        plt.axvline(0, c='k', ls='--')

    plot_name = alma_plot + alma_version + '/pattern_speeds_overview_alma'

    if mask_outside_bar:
        plot_name += '_bmask'
    else:
        plot_name += '_nobmask'

    if use_sim:
        plot_name += '_sim'

    plt.yticks(0.25 * np.array(range(len(alma_galaxies))), alma_galaxies)
    plt.xlabel(r'$\Omega_{p, \mathrm{TW}}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')

    plt.ylim([-0.125, 0.25 * len(alma_galaxies) - 0.125])

    if use_sim:
        x_lim = [-10, 10]
    else:
        x_lim = [0, 100]

    plt.xlim(x_lim)

    plt.tight_layout()

    plt.savefig(plot_name + '.png',
                bbox_inches='tight')
    plt.savefig(plot_name + '.pdf',
                bbox_inches='tight')

    plt.close()

    # Also create a KDE plot.

    plot_name = alma_plot + alma_version + '/pattern_speeds_kde_alma'

    if mask_outside_bar:
        plot_name += '_bmask'
    else:
        plot_name += '_nobmask'

    if use_sim:
        plot_name += '_sim'

    omega_bars = np.array(omega_bars)

    # Filter out any particularly weird pattern speeds (i.e. negative, very high)

    kde_mask = np.where((omega_bars < 100) & (omega_bars > 0))

    plt.figure(figsize=(8, 6))

    sns.kdeplot(omega_bars[kde_mask],
                bw='silverman', color='k', shade=True,
                clip=[0, 100])

    plt.ylabel('Probability Density')
    plt.xlabel(r'$\Omega_{p, \mathrm{TW}}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')

    if use_sim:
        x_lim = [-10, 10]
    else:
        x_lim = [0, 100]

    plt.xlim(x_lim)

    plt.savefig(plot_name + '.png',
                bbox_inches='tight')
    plt.savefig(plot_name + '.pdf',
                bbox_inches='tight')

    plt.close()

print('Complete!')
