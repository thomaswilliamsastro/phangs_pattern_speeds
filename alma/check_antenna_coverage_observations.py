# -*- coding: utf-8 -*-
"""
Check distributions of antenna coverage

@author: Tom Williams
"""

import os

from vars import phangs_folder

os.chdir(os.path.join(phangs_folder, 'alma'))

f = open('v34_download_log.txt', 'r')

combinations = {}

for line in f:
    if 'combination' in line:
        combination = line.split(':')[-1].strip()

        if combination in combinations.keys():
            combinations[combination] += 1
        else:
            combinations[combination] = 1

print(combinations)