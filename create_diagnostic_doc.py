# -*- coding: utf-8 -*-
"""
Put together a big PDF document for running through quality flagging

@author: Tom Williams
"""

import os
import time

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from vars import phangs_folder, alma_version, muse_version, plot_folder, alma_plot, muse_plot, alma_galaxies, \
    muse_galaxies, \
    hdu_types, star_masks, mask_outside_bars

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

if not os.path.exists(plot_folder + 'diagnostics'):
    os.mkdir(plot_folder + 'diagnostics')

start = time.time()

data_sources = ['alma', 'muse_mass', 'muse_ha']

# Assume the various parameters are the first in the list

hdu_type = hdu_types[0]
star_mask = star_masks[0]

for data_source in data_sources:

    print('Starting %s' % data_source)

    instrument = data_source.split('_')[0]

    if 'muse' in data_source:
        galaxies = muse_galaxies
        hdu_type = data_source.split('_')[-1]
    elif 'alma' in data_source:
        galaxies = alma_galaxies
    else:
        raise Warning('Weird data source!')

    if 'muse' in data_source:
        version = muse_version
    else:
        version = alma_version

    with PdfPages(plot_folder + 'diagnostics/' + data_source + '_' + version + '.pdf') as pdf:

        for galaxy in galaxies:

            galaxy = galaxy.strip()
            plt.figure(figsize=(8, 8))
            plt.suptitle(galaxy)

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

            try:
                bmask_im = plt.imread(bmask_file_name)
            except FileNotFoundError:
                plt.close()
                continue

            ax1 = plt.subplot(2, 2, 1)
            plt.imshow(bmask_im)
            plt.axis('off')
            plt.title('Beyond bar masked')

            # No bar masking

            ax2 = plt.subplot(2, 2, 2)

            nobmask_file_name = bmask_file_name.replace('bmask', 'nobmask')

            try:
                nobmask_im = plt.imread(nobmask_file_name)
                plt.imshow(nobmask_im)
                plt.title('No masking')
            except FileNotFoundError:
                plt.text(0.5, 0.5, 'Unmasked file not found',
                         ha='center', va='center',
                         transform=ax2.transAxes)

            plt.axis('off')

            # Slit lengths

            ax3 = plt.subplot(2, 2, 3)

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

            try:
                sl_im = plt.imread(sl_file_name)
                plt.imshow(sl_im)
                plt.title('Slit lengths')
            except FileNotFoundError:
                plt.text(0.5, 0.5, 'Slit length diagnostic not found',
                         ha='center', va='center',
                         transform=ax3.transAxes)

            plt.axis('off')

            # Slit lengths

            ax4 = plt.subplot(2, 2, 4)

            sw_file_name = sl_file_name.replace('slit_length', 'slit_width').replace('_sl', '_sw')

            try:
                sw_im = plt.imread(sw_file_name)
                plt.imshow(sw_im)
                plt.title('Slit widths')
            except FileNotFoundError:
                plt.text(0.5, 0.5, 'Slit width diagnostic not found',
                         ha='center', va='center',
                         transform=ax4.transAxes)

            plt.axis('off')

            plt.tight_layout()
            pdf.savefig(dpi=300)
            plt.close()

print('Complete! Took %.2fs' % (time.time() - start))
