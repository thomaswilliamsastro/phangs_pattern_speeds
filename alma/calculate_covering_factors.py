# -*- coding: utf-8 -*-
"""
Calculate ALMA covering factors to put into the diagnostic document

@author: Tom Williams
"""

import os
import warnings

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from vars import phangs_folder, alma_version, alma_galaxies, alma_output, galaxy_table

os.chdir(phangs_folder)

warnings.simplefilter('ignore')

galaxies = []
covering_factors = []

for galaxy in alma_galaxies:

    # From the table, pull out PA, inc, and R_e

    galaxy_row = galaxy_table[galaxy_table['name'] == galaxy]
    pa = galaxy_row['orient_posang'][0]
    inc = galaxy_row['orient_incl'][0]
    r_e = galaxy_row['size_reff'][0]
    ra = galaxy_row['orient_ra'][0]
    dec = galaxy_row['orient_dec'][0]

    galaxy = galaxy.upper()

    hdu_file_name = os.path.join('alma', alma_version, galaxy + '_mom0.fits')
    err_hdu_file_name = hdu_file_name.replace('mom0', 'emom0')

    try:
        hdu = fits.open(hdu_file_name)[0]
    except FileNotFoundError:
        continue
    err_hdu = fits.open(err_hdu_file_name)[0]

    w = WCS(hdu)
    pix_size = np.abs(hdu.header['CDELT1'] * 3600)
    x_cen, y_cen, _ = w.all_world2pix(ra, dec, 0, 1)

    grid_shape = hdu.data.shape
    grid_x = (np.arange(grid_shape[1]) - x_cen) * pix_size
    grid_y = (np.arange(grid_shape[0]) - y_cen) * pix_size

    x_coords, y_coords = np.meshgrid(grid_x, grid_y)

    angle = np.radians(pa)
    cos_a, sin_a = np.cos(angle), np.sin(angle)

    x_proj = x_coords * cos_a + y_coords * sin_a
    y_proj = - x_coords * sin_a + y_coords * cos_a

    # Account for inclination

    x_proj /= np.cos(np.radians(inc))

    x_proj = x_proj.reshape(grid_shape)
    y_proj = y_proj.reshape(grid_shape)
    r = np.sqrt(x_proj ** 2 + y_proj ** 2)

    ellipse = np.where(r <= 0.5 * r_e)

    # Find significant emission within r_e

    hdu.data[hdu.data <= err_hdu.data] = np.nan

    data_within_re = hdu.data[ellipse]

    all_data = len(data_within_re)
    significant_data = len(np.where(~np.isnan(data_within_re))[0])

    covering_factor = significant_data/all_data

    galaxies.append(galaxy)
    covering_factors.append(covering_factor)

cf_table = Table([galaxies, covering_factors],
                 names=['NAME', 'F_COV'])

cf_table.write(os.path.join(alma_output, 'covering_factors.fits'), overwrite=True)

print('Complete!')
