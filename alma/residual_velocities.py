# -*- coding: utf-8 -*-
"""
Look at the distribution of residual velocities from the Lang+ (2020) work

@author: Tom Williams
"""

import glob
import os
import warnings

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

import ps_functions
from vars import phangs_folder, galaxy_table, plot_folder

warnings.simplefilter('ignore')

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(os.path.join(phangs_folder, 'resids2020'))

resid_files = glob.glob('*.fits')

vel_cutoff = 50

all_resids = []
bar_resids = []
non_bar_resids = []

for resid_file in resid_files:
    hdu = fits.open(resid_file)[0]

    wcs = WCS(hdu)
    pix_size = np.abs(hdu.header['CDELT1'] * 3600)

    galaxy = resid_file.split('LOS')[0]

    phangs_row = galaxy_table[galaxy_table['name'] == galaxy]

    ra = phangs_row['orient_ra'][0]
    dec = phangs_row['orient_dec'][0]
    pa = phangs_row['orient_posang'][0]
    bar_r = phangs_row['morph_bar_r'][0]

    if np.isnan(bar_r):
        mask = np.zeros_like(hdu.data)
    else:
        x_cen, y_cen = wcs.all_world2pix(ra, dec, 1)

        grid_x = (np.arange(hdu.data.shape[1]) - x_cen) * pix_size
        grid_y = (np.arange(hdu.data.shape[0]) - y_cen) * pix_size

        x_coords, y_coords = np.meshgrid(grid_x, grid_y)

        bar_pix = ps_functions.bar_mask(x_coords, y_coords, pa, bar_r)

        mask = np.ones_like(hdu.data)
        mask[bar_pix] = 0

    # Pull out all of the residual velocities
    data_flat = hdu.data.flatten()
    data_flat = data_flat[~np.isnan(data_flat)]
    data_flat = data_flat[np.abs(data_flat) <= vel_cutoff]
    all_resids.extend(data_flat)

    bar_flat = hdu.data[mask == 1].flatten()
    bar_flat = bar_flat[~np.isnan(bar_flat)]
    bar_flat = bar_flat[np.abs(bar_flat) <= vel_cutoff]
    bar_resids.extend(bar_flat)

    non_bar_flat = hdu.data[mask == 0].flatten()
    non_bar_flat = non_bar_flat[~np.isnan(non_bar_flat)]
    non_bar_flat = non_bar_flat[np.abs(non_bar_flat) <= vel_cutoff]
    non_bar_resids.extend(non_bar_flat)

os.chdir('../')
plot_name = os.path.join(plot_folder, 'residual_hist')

plt.figure(figsize=(5, 3))
_, bins, _ = plt.hist(all_resids, bins=20, color='k', histtype='step', label='All')
plt.hist(non_bar_resids, bins=bins, color='b', histtype='step', hatch='/', label='Non-bar')
plt.hist(bar_resids, bins=bins, color='r', histtype='step', hatch='\\', label='Bar')

plt.xlabel(r'$v_\mathrm{residual}$ (km s$^{-1}$)')
plt.ylabel(r'$N_\mathrm{pix}$')

plt.legend(loc='upper left', frameon=False, framealpha=0.5)

plt.tight_layout()
plt.savefig(plot_name + '.png', bbox_inches='tight')
plt.savefig(plot_name + '.pdf', bbox_inches='tight')
plt.close()

print('Complete!')
