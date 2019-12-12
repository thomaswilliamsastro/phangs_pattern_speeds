# -*- coding: utf-8 -*-
"""
Download intensity/velocity maps and errors for all available ALMA data at highest possible resolution

@author: Tom Williams
"""

import getpass
import os

import pexpect
from astropy.io import fits
from astropy.table import Table
from tqdm import tqdm

from alma.folders import phangs_folder

# TODO: Update to new ALMA data release.

os.chdir(phangs_folder)

overwrite_data = False

username = 'williams'
password = getpass.getpass()
server = 'astro-node4'
data_folder = 'data/beegfs/astro-storage/groups/schinnerer/PHANGS/ALMA/Data_Release/' \
              'PHANGS-ALMA-v33-by-Adam/delivery/broad_maps'

# Open up the table for latest data release (as of 20191211)

galaxy_table = fits.open('documents/phangs_sample_table_v1p2.fits')
galaxy_table = Table(galaxy_table[1].data)

galaxies = galaxy_table['NAME'][galaxy_table['ALMA'] == 1]

# For each galaxy, download the ALMA data as appropriate. Use broad maps.

for galaxy in tqdm(galaxies):

    # Construct the command to download ALMA data from the server. We want mom0, mom1 (and associated errors)

    galaxy = galaxy.strip()
    gal_low = galaxy.lower()

    # Some pesky galaxies are saved under aliases, not the NAME column in the table. Force these in manually.

    try:
        gal_low = {'ESO097-013': 'circinus',
                   }[galaxy]
    except KeyError:
        pass

    print('Downloading ALMA data for ' + galaxy)

    # Setups in order of preference (most combinations, highest resolutions first).

    antenna_setups = ['12m+7m+tp', '12m+7m', '7m+tp', '7m']

    for antenna_setup in antenna_setups:

        file_name = gal_low + '_' + antenna_setup + '_co21_broad_mom0.fits'

        command = ('ssh '
                   + username + '@' + server + ' test -e /' + data_folder + '/' + file_name
                   + '; echo result: $?'
                   )

        var_child = pexpect.spawn(command)
        var_child.expect(["password:", pexpect.EOF])
        var_child.sendline(password)

        result = var_child.expect(['result: 0', 'result: 1'], timeout=600)

        if result == 0:
            break

    extensions = ['mom0', 'emom0', 'mom1', 'emom1']

    for extension in extensions:

        if not os.path.exists('alma/' + galaxy + '_' + extension + '.fits') or overwrite_data:
            file_name = gal_low + '_' + antenna_setup + '_co21_broad_' + extension + '.fits'

            command = ('scp '
                       + username + '@' + server + ':/' + data_folder + '/' + file_name
                       + ' alma/' + galaxy + '_' + extension + '.fits')

            # SCP requires a password (input earlier), so use pexpect to pass that along. Use a 5 min timeout for large
            # files

            var_child = pexpect.spawn(command)
            var_child.expect(["password:", pexpect.EOF])
            var_child.sendline(password)
            var_child.expect(pexpect.EOF, timeout=600)

print('Complete!')
