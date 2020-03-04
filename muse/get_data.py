# -*- coding: utf-8 -*-
"""
Download data for each galaxy that has an environmental mask

@author: Tom Williams
"""

import getpass
import os

import pexpect

from vars import phangs_folder, muse_galaxies

os.chdir(phangs_folder)

if not os.path.exists('muse'):
    os.mkdir('muse')

overwrite_data = False

username = 'williams'
password = getpass.getpass()
server = 'astro-node4'
data_folder = '/data/beegfs/astro-storage/groups/schinnerer/pessa/DR1/data'

# For each galaxy, download the MUSE data as appropriate.

# TODO: When DR2 appears, update. Refactor this code to version the files into the right folder.

for galaxy in muse_galaxies:

    # Construct the command to download MUSE data from the server.

    if not os.path.exists('muse/' + galaxy + '_MAPS.fits') or overwrite_data:
        print('Downloading MUSE data for ' + galaxy)

        file_name = galaxy + '/' + galaxy + '_v1' + '/' + galaxy + '_MAPS.fits'

        command = ('scp '
                   + username + '@' + server + ':/' + data_folder + '/' + file_name
                   + ' ./muse/')

        # SCP requires a password (input earlier), so use pexpect to pass that along. Because some of this data is
        # quite Lorge, use a 1hr timeout.

        var_child = pexpect.spawn(command)
        var_child.expect(["password:", pexpect.EOF])
        var_child.sendline(password)
        var_child.expect(pexpect.EOF, timeout=3600)

print('Complete!')
