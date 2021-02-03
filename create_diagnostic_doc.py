# -*- coding: utf-8 -*-
"""
Put together a big PDF document for running through quality flagging

@author: Tom Williams
"""

import os
import time

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages

from vars import phangs_folder, alma_version, muse_version, plot_folder, alma_plot, muse_plot, alma_galaxies, \
    hdu_types, star_masks, mask_outside_bars, alma_output

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

if not os.path.exists(plot_folder + 'diagnostics'):
    os.mkdir(plot_folder + 'diagnostics')

start = time.time()

data_sources = ['alma', 'muse_mass', 'muse_ha']

f_cov_table = Table.read(alma_output + 'covering_factors.fits')

# Assume the various parameters are the first in the list

hdu_type = hdu_types[0]
star_mask = star_masks[0]

with PdfPages(plot_folder + 'diagnostics/diagnostics.pdf') as pdf:

    for galaxy in alma_galaxies:

        alma_found = False
        muse_found = False

        galaxy = galaxy.upper()

        # print('Beginning %s' % galaxy)

        galaxy_dict = {}

        for data_source in data_sources:

            instrument = data_source.split('_')[0]

            if 'muse' in data_source:
                hdu_type = data_source.split('_')[-1]

            if 'muse' in data_source:
                version = muse_version
            else:
                version = alma_version

            # Beyond bar masking

            if 'muse' in data_source:
                bmask_file_name = muse_plot + muse_version + '/'
            else:
                bmask_file_name = alma_plot + alma_version + '/'

            if 'muse' in data_source:
                bmask_file_name += hdu_type + '_'
                if star_mask:
                    bmask_file_name += 'smask_'
                else:
                    bmask_file_name += 'nosmask_'

            if mask_outside_bars:
                bmask_file_name += 'bmask/'
            else:
                bmask_file_name += 'nobmask/'

            bmask_file_name += galaxy + '_'

            if 'muse' in data_source:
                bmask_file_name += hdu_type + '_'
                if star_mask:
                    bmask_file_name += 'smask_'
                else:
                    bmask_file_name += 'nosmask_'

            if mask_outside_bars:
                bmask_file_name += 'bmask_'
            else:
                bmask_file_name += 'nobmask_'

            bmask_file_name += instrument + '.png'
            nobmask_file_name = bmask_file_name.replace('bmask', 'nobmask')

            if 'muse' in data_source:
                sl_file_name = muse_plot + muse_version + '/'
            else:
                sl_file_name = alma_plot + alma_version + '/'

            sl_file_name += 'slit_length/'

            sl_file_name += galaxy + '_'

            if 'muse' in data_source:
                sl_file_name += hdu_type + '_'
                if star_mask:
                    sl_file_name += 'smask_'
                else:
                    sl_file_name += 'nosmask_'

            if mask_outside_bars:
                sl_file_name += 'bmask_'
            else:
                sl_file_name += 'nobmask_'

            sl_file_name += instrument + '_sl_comparison.png'
            sw_file_name = sl_file_name.replace('slit_length', 'slit_width').replace('_sl', '_sw')

            try:
                bmask_im = plt.imread(bmask_file_name)
                galaxy_dict[data_source + '_bmask'] = bmask_im
                galaxy_dict[data_source + '_nobmask'] = plt.imread(nobmask_file_name)
                galaxy_dict[data_source + '_sl'] = plt.imread(sl_file_name)
                galaxy_dict[data_source + '_sw'] = plt.imread(sw_file_name)

                if data_source == 'alma':
                    galaxy_dict[data_source + '_strict'] = plt.imread(bmask_file_name.replace('_bmask',
                                                                                              '_bmask_strict'))
                    galaxy_dict[data_source + '_lowres'] = plt.imread(bmask_file_name.replace('_bmask',
                                                                                              '_bmask_lowres'))

                if 'muse' in data_source:
                    muse_found = True
                if 'alma' in data_source:
                    f_cov = f_cov_table[f_cov_table['NAME'] == galaxy]['F_COV'][0]
                    alma_found = True

                    f = open('alma/v34_download_log.txt', 'r')

                    for line in f:
                        if galaxy + ' antenna combination' in line:
                            antenna_combination = line.split(':')[-1].strip()

                    f.close()

                    galaxy_dict['alma_config'] = antenna_combination
            except FileNotFoundError:
                pass

        # print(alma_found, muse_found)
        if not alma_found:
            f_cov = 0

        if not alma_found and not muse_found:
            print('%s not found' % galaxy)
            continue

        n_panels = np.sum([alma_found, muse_found, muse_found])

        fig_width = n_panels * 8

        n_cols = n_panels * 2
        n_rows = 3

        plt.figure(figsize=(fig_width, 14))
        plt.suptitle('%s (%s, f_cov=%.2f)' % (galaxy, antenna_combination, f_cov))

        n_subplot = 1

        fancy_titles = {'alma': 'ALMA',
                        'muse_ha': r'MUSE H$\alpha$',
                        'muse_mass': r'MUSE M$_\ast$'}[data_source]

        final_data_sources = []

        if alma_found:
            final_data_sources.extend(['alma'])

        if muse_found:
            final_data_sources.extend(['muse_mass', 'muse_ha'])

        for data_source in final_data_sources:
            fancy_title = {'alma': 'ALMA',
                           'muse_ha': r'MUSE H$\alpha$',
                           'muse_mass': r'MUSE M$_\ast$'}[data_source]

            plt.subplot(n_rows, n_cols, n_subplot)
            plt.imshow(galaxy_dict[data_source + '_bmask'])
            plt.axis('off')
            plt.title('%s, beyond bar masked' % fancy_title)

            n_subplot += 1

            plt.subplot(n_rows, n_cols, n_subplot)
            plt.imshow(galaxy_dict[data_source + '_nobmask'])
            plt.axis('off')
            plt.title('%s, no masking' % fancy_title)

            n_subplot += 1

        for data_source in final_data_sources:
            plt.subplot(n_rows, n_cols, n_subplot)
            plt.imshow(galaxy_dict[data_source + '_sl'])
            plt.title('Slit lengths')
            plt.axis('off')

            n_subplot += 1

            plt.subplot(n_rows, n_cols, n_subplot)
            plt.imshow(galaxy_dict[data_source + '_sw'])
            plt.title('Slit widths')
            plt.axis('off')

            n_subplot += 1

        # Finally, add in the strict map, and the low resolution map for the ALMA data

        if alma_found:

            plt.subplot(n_rows, n_cols, n_subplot)

            plt.imshow(galaxy_dict['alma_strict'])
            plt.title('Strict')
            plt.axis('off')

            n_subplot += 1

            plt.subplot(n_rows, n_cols, n_subplot)

            plt.imshow(galaxy_dict['alma_lowres'])
            plt.title('Low-res')
            plt.axis('off')

        plt.tight_layout()
        # plt.show()
        pdf.savefig(dpi=300)
        plt.close()

print('Complete! Took %.2fs' % (time.time() - start))
