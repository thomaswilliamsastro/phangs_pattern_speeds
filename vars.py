# -*- coding: utf-8 -*-
"""
Frequently used variables in the code

@author: Tom Williams
"""

import socket

import numpy as np
from astropy.io import fits
from astropy.table import Table

if socket.gethostname() in ['astro-node4']:
    phangs_folder = '/data/beegfs/astro-storage/groups/schinnerer/williams/phangs'
else:
    phangs_folder = '/Users/williams/Documents/phangs'

output_folder = 'pattern_speeds_output/'
alma_output = output_folder + 'alma/'
muse_output = output_folder + 'muse/'
resonance_output = output_folder + 'resonance/'

plot_folder = 'plots/pattern_speeds/'
alma_plot = plot_folder + 'alma/'
muse_plot = plot_folder + 'muse/'
resonance_plot = plot_folder + 'resonance/'
radial_plot = plot_folder + 'radial_profiles/'

megatable_folder = 'mega_tables/'

# Version numbering for things

table_version = 'v1p5'
alma_version = 'v34'
muse_version = 'DR2.0'
pattern_speed_version = 'v1p0'

# Galaxies to fit

galaxy_table = fits.open(phangs_folder + '/documents/phangs_sample_table_' + table_version + '.fits')
galaxy_table = Table(galaxy_table[1].data)

# Pull out galaxies that have ALMA data

alma_galaxies = galaxy_table['name'][galaxy_table['survey_alma_status'] == 'released']
# alma_galaxies = galaxy_table['NAME'][galaxy_table['HAS_ALMA'] == 1]
alma_galaxies = sorted(alma_galaxies)

# TODO: Pull out galaxies that have MUSE data

# muse_galaxies = galaxy_table['name'][galaxy_table['survey_muse_status'] == 'released']

muse_galaxies = ['ic5332', 'ngc0628', 'ngc1087', 'ngc1300', 'ngc1365', 'ngc1385', 'ngc1433', 'ngc1512', 'ngc1566',
                 'ngc1672', 'ngc2835', 'ngc3351', 'ngc3627', 'ngc4254', 'ngc4303', 'ngc4321', 'ngc4535', 'ngc5068',
                 'ngc7496']
muse_galaxies = sorted(muse_galaxies)

# Various options we can play around with for testing

overwrite_pafit = False
overwrite_bootstraps = False

mask_outside_bars = [True]
star_masks = [True]
slit_widths = [1]
slit_lengths = [0]
hdu_types = ['mass', 'ha']  # ['mass', 'spitzer_mass', 'flux', 'ha', 'toy_sim', 'toy_sim_cov']
hdu_types = ['mass']

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

# To stop large number of plots, turn off plotting if we're varying slit lengths/slit widths or masking with GMCs

use_gmc_mask = False

if use_gmc_mask:
    gmc_removal_factors = np.arange(0.1, 1.0, 0.1)
    alma_galaxies = ['ngc0300', 'ngc1512', 'ngc2090', 'ngc3351']
else:
    gmc_removal_factors = [0]

if len(slit_lengths) > 1 or len(slit_widths) > 1:
    plot = False
else:
    plot = True
