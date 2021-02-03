# -*- coding: utf-8 -*-
"""
Do a nice plot of resonances

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

from vars import phangs_folder, output_folder, muse_version, pattern_speed_version, galaxy_table, \
    plot_folder


def get_data_and_resonances(galaxy_table, pattern_speeds_table, galaxy):
    table_row = galaxy_table[galaxy_table['name'] == galaxy]
    inc = table_row['orient_incl'][0]
    pa = table_row['orient_posang'][0]
    dist = table_row['dist'][0]
    ra = table_row['orient_ra'][0]
    dec = table_row['orient_dec'][0]

    galaxy = galaxy.upper()

    # Pull in corotation radius

    pattern_speeds_row = pattern_speeds_table[pattern_speeds_table['GALAXY'] == galaxy]
    r_cr = pattern_speeds_row['R_CR_MUSE_MASS_1'][0]
    r_cr_err = pattern_speeds_row['R_CR_MUSE_MASS_1_ERR'][0]
    r_ilr = pattern_speeds_row['R_ILR_MUSE_MASS_1'][0]
    r_ilr_err = pattern_speeds_row['R_ILR_MUSE_MASS_1_ERR'][0]
    r_olr = pattern_speeds_row['R_OLR_MUSE_MASS_1'][0]
    r_olr_err = pattern_speeds_row['R_OLR_MUSE_MASS_1_ERR'][0]

    r_cr_upper = r_cr + r_cr_err
    r_cr_lower = r_cr - r_cr_err

    r_ilr_upper = r_ilr + r_ilr_err
    r_ilr_lower = r_ilr - r_ilr_err

    r_olr_upper = r_olr + r_olr_err
    r_olr_lower = r_olr - r_olr_err

    hdu_file_name = os.path.join('muse', muse_version, galaxy + '_MAPS.fits')

    hdu = fits.open(hdu_file_name)['FLUX']
    data = np.log10(hdu.data)
    w = WCS(hdu)
    pix_size = np.abs(hdu.header['CD1_1']) * 3600

    x_cen, y_cen = w.all_world2pix(ra, dec, 1)

    xi, yi = np.meshgrid((np.arange(data.shape[1]) - x_cen),
                         (np.arange(data.shape[0]) - y_cen))

    # Convert these positions to physical positions (kpc), accounting for inclination and rotation

    xi *= pix_size / 3600 * np.pi / 180 * dist * 1e3
    yi *= pix_size / 3600 * np.pi / 180 * dist * 1e3

    # Create projected radius

    angle = np.radians(pa)
    cos_a, sin_a = np.cos(angle), np.sin(angle)

    x_proj = xi * cos_a + yi * sin_a
    y_proj = - xi * sin_a + yi * cos_a

    # Account for inclination

    x_proj /= np.cos(np.radians(inc))

    r = np.sqrt(x_proj ** 2 + y_proj ** 2)

    resonances = np.zeros_like(r)
    resonances[resonances == 0] = np.nan

    resonances[(r <= r_cr_upper) & (r >= r_cr_lower)] = 1
    resonances[(r <= r_ilr_upper) & (r >= r_ilr_lower)] = 0
    resonances[(r <= r_olr_upper) & (r >= r_olr_lower)] = 0

    return data, resonances, w, dist


matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

pattern_speeds_table = Table.read(output_folder + 'pattern_speed_table_' + pattern_speed_version + '.fits')

good_galaxies = pattern_speeds_table[(pattern_speeds_table['OM_P_MUSE_MASS_QUAL'] == 1) |
                                     (pattern_speeds_table['OM_P_MUSE_MASS_QUAL'] == 2)]['GALAXY']

good_galaxies = ['ngc3351']

for galaxy in good_galaxies:
    data, resonances, w, dist = get_data_and_resonances(galaxy_table, pattern_speeds_table, galaxy.lower())

    kpc_scale = 1 / (0.2 / 3600 * np.pi / 180 * dist * 1e3)

    vmin = np.nanpercentile(data, 0.5)
    vmax = np.nanpercentile(data, 99.5)

    plot_name = os.path.join(plot_folder, 'resonance', galaxy + '_resonance_map')

    plt.figure(figsize=(4, 4))
    ax = plt.subplot(projection=w)

    plt.imshow(data, origin='lower', cmap='viridis', vmin=vmin, vmax=vmax)
    plt.imshow(resonances, origin='lower', cmap='rainbow', alpha=0.5, vmin=0, vmax=1)

    plt.xlabel('RA (J2000)')
    plt.ylabel('Dec (J2000)')

    scalebar = AnchoredSizeBar(ax.transData,
                               kpc_scale, '1kpc', 3,
                               pad=0.3,
                               borderpad=0.3,
                               color='black',
                               frameon=True)

    ax.add_artist(scalebar)

    plt.tight_layout()

    # plt.show()
    plt.savefig(plot_name + '.pdf', bbox_inches='tight')
    plt.savefig(plot_name + '.png', bbox_inches='tight')
    plt.close()

print('Complete!')
