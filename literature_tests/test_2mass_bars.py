# -*- coding: utf-8 -*-
"""
Test correlations with 2MASS bar sample

@author: Tom Williams
"""

import os

import matplotlib
import matplotlib.pyplot as plt

from vars import phangs_folder

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder + '/lit_tables')

# Read in the inclinations and galaxy names

gal_info = {}

f = open('2mass_tab_1.txt', 'r')

for i in range(30):
    f.readline()

for line in f:
    galaxy = line[5:13].strip()
    inclination = float(line[62:64])

    gal_info[galaxy] = [inclination]

# Read in the deprojected bar lengths

f = open('2mass_tab_2.txt', 'r')

for i in range(25):
    f.readline()

for line in f:

    galaxy_found = False

    for galaxy in gal_info.keys():

        if line[23:31].strip() == galaxy:
            gal_info[galaxy].append(float(line[76:80]))

f.close()

inclinations = []
bar_r_deproj = []

for key in gal_info.keys():

    if len(gal_info[key]) == 2:
        inclinations.append(gal_info[key][0])
        bar_r_deproj.append(gal_info[key][1])

plt.figure(figsize=(8, 6))

plt.scatter(inclinations, bar_r_deproj, c='k')

plt.xlabel(r'$i$ (deg)')
plt.ylabel(r'$R_\mathrm{bar}$ (kpc)')

plt.savefig('2_mass_r_deproj.png',
            bbox_inches='tight')

print('Complete!')
