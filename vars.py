# -*- coding: utf-8 -*-
"""
Frequently used variables in the code

@author: Tom Williams
"""

from astropy.io import fits
from astropy.table import Table

phangs_folder = '/Users/williams/Documents/phangs'

output_folder = 'pattern_speeds_output/'
alma_output = output_folder + 'alma/'
muse_output = output_folder + 'muse/'
corot_output = output_folder + 'corot/'

plot_folder = 'plots/pattern_speeds/'
alma_plot = plot_folder + 'alma/'
muse_plot = plot_folder + 'muse/'
corot_plot = plot_folder + 'corot/'
radial_plot = plot_folder + 'radial_profiles/'

megatable_folder = 'mega_tables/'

# Version numbering for things

table_version = 'v1p4'
alma_version = 'v34'
muse_version = 'DR1.0'
corot_version = 'v0p1'

# Galaxies to fit

galaxy_table = fits.open(phangs_folder + '/documents/phangs_sample_table_' + table_version + '.fits')
galaxy_table = Table(galaxy_table[1].data)

# Pull out galaxies that have ALMA data

alma_galaxies = galaxy_table['NAME'][galaxy_table['HAS_ALMA'] == 1]
alma_galaxies = sorted(alma_galaxies)

muse_galaxies = ['NGC1087', 'NGC1512', 'NGC1672', 'NGC3351', 'NGC4254', 'NGC5068',
                 'IC5332', 'NGC1365', 'NGC1566', 'NGC2835', 'NGC3627', 'NGC4535', 'NGC628']
muse_galaxies = sorted(muse_galaxies)

# TODO: Pull out galaxies that have MUSE data

# muse_galaxies = galaxy_table['NAME'][galaxy_table['HAS_MUSE'] == 1]
# muse_galaxies = sorted(muse_galaxies)

# Various options we can play around with for testing

overwrite_pafit = False
overwrite_bootstraps = False

mask_outside_bars = [True]
star_masks = [True]
slit_widths = [1]
slit_lengths = [0]
hdu_types = ['mass']  # ['mass', 'flux', 'ha', 'toy_sim', 'toy_sim_cov']

# Number of pattern speeds to fit (really only 1 or 2)

n_pattern_speeds = 1

# If n_pattern_speeds > 1, we shouldn't mask beyond the bar

if n_pattern_speeds > 1:
    mask_outside_bars = [False]

# Switch for simulation testing for FOV issues

use_sim = False

# Slit lengths

# start = 10
# stop = 150
# step = 10
#
# slit_lengths = np.arange(start, stop + step, step)

# Slit widths

# start = 1
# stop = 10
# step = 0.5
#
# slit_widths = np.arange(start, stop + step, step)

# To stop large number of plots, turn off plotting if we're varying slit lengths/slit widths

if len(slit_lengths) > 1 or len(slit_widths) > 1:
    plot = False
else:
    plot = True
