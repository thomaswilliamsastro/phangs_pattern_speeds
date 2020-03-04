# -*- coding: utf-8 -*-
"""
Look at radial profiles from the Mega-Tables

@author: Tom Williams
"""

import os
import time

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table

from vars import phangs_folder, corot_version, output_folder, radial_plot, megatable_folder, alma_galaxies

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

start = time.time()

os.chdir(phangs_folder)

if not os.path.exists(radial_plot):
    os.mkdir(radial_plot)

# Read in the pattern speed table

pattern_speed_table = fits.open(output_folder + 'pattern_speed_table_' + corot_version + '.fits')
pattern_speed_table = Table(pattern_speed_table[1].data)

alma_corot_colnames = []
muse_corot_colnames = []

for colname in pattern_speed_table.colnames:

    if 'R_CR' in colname and 'ERR' not in colname:

        if 'ALMA' in colname:
            alma_corot_colnames.append(colname)

        elif 'MUSE' in colname:
            muse_corot_colnames.append(colname)

# alma_galaxies = ['NGC1559']

quantity_colnames = ['Sigma_mol']

for galaxy in alma_galaxies:

    has_alma = False
    has_muse = False

    galaxy = galaxy.strip()

    print('Beginning %s' % galaxy)

    # For each galaxy, pull out the calculated corotation radii (both MUSE and ALMA)

    row = pattern_speed_table[pattern_speed_table['GALAXY'] == galaxy]

    if len(row) == 0:
        print('No pattern speeds calculated: skipping')
        continue

    alma_corot_r = np.array([row[colname][0] for colname in alma_corot_colnames])
    alma_corot_r_err = np.array([row[colname + '_ERR'][0] for colname in alma_corot_colnames])

    if np.any(np.isnan(alma_corot_r) == False):
        has_alma = True

    muse_corot_r = np.array([row[colname][0] for colname in muse_corot_colnames])
    muse_corot_r_err = np.array([row[colname + '_ERR'][0] for colname in muse_corot_colnames])

    if np.any(np.isnan(muse_corot_r) == False):
        has_muse = True

    if np.all(np.isnan(alma_corot_r)) and np.all(np.isnan(muse_corot_r)):
        print('No corotations found: skipping')
        continue

    # Read in the radial profile megatable for the galaxy

    galaxy_megatable = fits.open(megatable_folder + galaxy + '_radial_profile_200pc_phys.fits')
    galaxy_megatable = Table(galaxy_megatable[1].data)

    radius = (galaxy_megatable['r_gal_min'] + galaxy_megatable['r_gal_max']) / 2

    for quantity_colname in quantity_colnames:
        x_label = r'$R$ (kpc)'

        try:
            y_label = {'Sigma_mol': r'$\Sigma_\mathrm{mol}$ (M$_\odot$ pc$^{-2}$)',
                       }[quantity_colname]
        except KeyError:
            y_label = quantity_colname
            print('No label defined for %s -- tidy this up!' % quantity_colname)

        quantity = galaxy_megatable[quantity_colname]

        if np.all(np.isnan(quantity)):
            print('No radial profile found: skipping')
            continue

        plt.figure(figsize=(8, 6))

        plt.scatter(radius, quantity,
                    c='k')

        y_min, y_max = plt.ylim()

        if has_alma:
            plt.vlines(alma_corot_r, y_min, y_max, colors='b', label='ALMA')
            plt.vlines(alma_corot_r + alma_corot_r_err, y_min, y_max, colors='b', linestyles='dashed')
            plt.vlines(alma_corot_r - alma_corot_r_err, y_min, y_max, colors='b', linestyles='dashed')

        if has_muse:
            plt.vlines(muse_corot_r, y_min, y_max, colors='r', label='MUSE')
            plt.vlines(muse_corot_r + muse_corot_r_err, y_min, y_max, colors='r', linestyles='dashed')
            plt.vlines(muse_corot_r - muse_corot_r_err, y_min, y_max, colors='r', linestyles='dashed')

        plt.xlabel(x_label)
        plt.ylabel(y_label)

        plt.ylim([y_min, y_max])

        plt.legend(loc='upper right', frameon=False)

        if not os.path.exists(radial_plot + quantity_colname):
            os.mkdir(radial_plot + quantity_colname)

        plt.savefig(radial_plot + quantity_colname + '/' + galaxy + '_' + quantity_colname + '.png',
                    bbox_inches='tight')
        plt.savefig(radial_plot + quantity_colname + '/' + galaxy + '_' + quantity_colname + '.pdf',
                    bbox_inches='tight')

        plt.close()

print('Complete! Took %.2fs' % (time.time()-start))
