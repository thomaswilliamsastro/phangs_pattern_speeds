# -*- coding: utf-8 -*-
"""
Test effect of masking clouds on recovered pattern speeds

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

from vars import phangs_folder, output_folder, plot_folder, gmc_removal_factors, alma_galaxies, corot_version

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

corot_table = Table.read(output_folder + 'pattern_speed_table_' + corot_version + '.fits')

for galaxy in alma_galaxies:

    galaxy = galaxy.upper()

    corot_row = corot_table[corot_table['GALAXY'] == galaxy]
    true_om_p = corot_row['OM_P_ALMA'][0]
    true_om_p_err_up = corot_row['OM_P_ALMA_ERR_UP'][0]
    true_om_p_err_down = corot_row['OM_P_ALMA_ERR_DOWN'][0]

    omega_ps = []
    omega_ps_err_up = []
    omega_ps_err_down = []

    plot_filename = plot_folder + 'gmc/'

    if not os.path.exists(plot_filename):
        os.makedirs(plot_filename)

    plot_filename += galaxy + '_gmc_comparison'

    for gmc_removal_factor in gmc_removal_factors:

        try:

            om_p, om_p_err_up, om_p_err_down = np.loadtxt('%salma/v34/gmc_mask/%s_bmask_gmc_%.1f_pattern_speed_alma.txt'
                                                          % (output_folder, galaxy, gmc_removal_factor),
                                                          unpack=True)

        except OSError:
            om_p, om_p_err_up, om_p_err_down = np.nan, np.nan, np.nan

        omega_ps.append(om_p)
        omega_ps_err_up.append(om_p_err_up)
        omega_ps_err_down.append(om_p_err_down)

    omega_ps = np.array(omega_ps)
    omega_ps_err_up = np.array(omega_ps_err_up)
    omega_ps_err_down = np.array(omega_ps_err_down)

    plt.figure(figsize=(4, 3))

    plt.errorbar(gmc_removal_factors * 100, omega_ps, yerr=[omega_ps_err_down, omega_ps_err_up],
                 c='k', ls='none', marker='o')

    plt.axhline(true_om_p, c='k', ls='-')
    plt.axhline(true_om_p + true_om_p_err_up, c='k', ls='--')
    plt.axhline(true_om_p - true_om_p_err_down, c='k', ls='--')

    plt.ylabel(r'$\Omega_{p}\, (\mathrm{km\,s}^{-1}\,\mathrm{kpc}^{-1})$')
    plt.xlabel('GMC Removal (%)')

    plt.tight_layout()

    # plt.title(galaxy)

    # plt.show()

    plt.savefig(plot_filename + '.png', bbox_inches='tight')
    plt.savefig(plot_filename + '.pdf', bbox_inches='tight')

    plt.close()

print('Complete!')
