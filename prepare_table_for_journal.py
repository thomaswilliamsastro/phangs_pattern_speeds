# -*- coding: utf-8 -*-
"""
Trim the fits table down to 3 DP for the resonance columns

@author: Tom Williams
"""

import os

from astropy.table import Table

from vars import phangs_folder, output_folder, pattern_speed_version

os.chdir(os.path.join(phangs_folder, output_folder))

table = Table.read('pattern_speed_table_%s.fits' % pattern_speed_version)

names_to_trim = ['OM_P', 'R_CR', 'R_ILR', 'R_OLR']
names_to_skip = ['REF', 'QUAL']
colnames = table.colnames

for colname in colnames:

    if any(col in colname for col in names_to_skip):
        continue

    if any(col in colname for col in names_to_trim):
        table[colname] = table[colname].round(3)

table.write('pattern_speed_table_%s_edit.fits' % pattern_speed_version, overwrite=True)
print('complete!')
