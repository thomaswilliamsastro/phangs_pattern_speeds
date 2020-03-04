# -*- coding: utf-8 -*-
"""
Download intensity/velocity maps and errors for all available ALMA data at highest possible resolution

@author: Tom Williams
"""

import getpass
import os
import time

import pexpect

from vars import phangs_folder, alma_galaxies, alma_version

start = time.time()

os.chdir(phangs_folder)

if not os.path.exists('alma/'):
    os.mkdir('alma/')
if not os.path.exists('alma/' + alma_version):
    os.mkdir('alma/' + alma_version)

overwrite_data = False

username = 'williams'
password = getpass.getpass()
server = 'astro-node4'
server_folder = 'data/beegfs/astro-storage/groups/schinnerer/PHANGS/ALMA/Data_Release/'
data_folder = 'PHANGS-ALMA-' + alma_version + '-by-Adam/broad_maps'

# Set up a log so we can refer back to the particular antenna setup

log = open('alma/' + alma_version + '_download_log.txt', 'w+')

# For each galaxy, download the ALMA data as appropriate. Use broad maps.

for i, galaxy in enumerate(alma_galaxies):

    # Construct the command to download ALMA data from the server. We want mom0, mom1 (and associated errors)

    galaxy = galaxy.strip()
    gal_low = galaxy.lower()

    # Some pesky galaxies are saved under aliases, not the NAME column in the table. Force these in manually.

    try:
        gal_low = {'ESO097-013': 'circinus',
                   }[galaxy]
    except KeyError:
        pass

    log.write('%d/%d: %s\n' % (i+1, len(alma_galaxies), galaxy))
    print('%d/%d: %s' % (i+1, len(alma_galaxies), galaxy))

    # Setups in order of preference (most combinations, highest resolutions first).

    antenna_setups = ['12m+7m+tp', '12m+7m', '7m+tp', '7m']
    setup_found = False

    for antenna_setup in antenna_setups:

        file_name = gal_low + '_' + antenna_setup + '_co21_broad_mom0.fits'

        command = ('ssh '
                   + username + '@' + server + ' test -e /' + server_folder + data_folder + '/' + file_name
                   + '; echo result: $?'
                   )

        var_child = pexpect.spawn(command)
        var_child.expect(["password:", pexpect.EOF])
        var_child.sendline(password)

        result = var_child.expect(['result: 0', 'result: 1'], timeout=600)

        if result == 0:
            setup_found = True
            break

    if not setup_found:
        log.write('No suitable combination found for %s: skipping\n' % galaxy)
    else:
        log.write('%s antenna combination: %s\n' % (galaxy, antenna_setup))

    extensions = ['mom0', 'emom0', 'mom1', 'emom1']

    for extension in extensions:

        if not os.path.exists('alma/' + alma_version + '/' + galaxy + '_' + extension + '.fits') or overwrite_data:
            file_name = gal_low + '_' + antenna_setup + '_co21_broad_' + extension + '.fits'

            command = ('scp '
                       + username + '@' + server + ':/' + server_folder + data_folder + '/' + file_name
                       + ' alma/' + alma_version + '/' + galaxy + '_' + extension + '.fits')

            # SCP requires a password (input earlier), so use pexpect to pass that along. Use a 5 min timeout for large
            # files

            var_child = pexpect.spawn(command)
            var_child.expect(["password:", pexpect.EOF])
            var_child.sendline(password)
            var_child.expect(pexpect.EOF, timeout=600)

log.write('Complete! Took %.2fm\n' % ((time.time()-start)/60))

log.close()

print('Complete!')
