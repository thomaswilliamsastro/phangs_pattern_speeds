# -*- coding: utf-8 -*-
"""
Calculate the corotation radius for the PHANGS MUSE galaxies

@author: Tom Williams
"""

import os
import time

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

from vars import phangs_folder, alma_version, muse_version, alma_output, muse_output, alma_galaxies, galaxy_table, \
    star_masks, mask_outside_bars, corot_output, corot_plot


def find_corotation_radii(om, om_up, om_down, pattern_speed, pattern_speed_up, pattern_speed_down,
                          sigma=1):
    om_corot = []
    r_corot = []
    corotation = {}

    idx = 0

    for i, _ in enumerate(om):

        diff_up_up = (om[i] + sigma * om_up[i]) - (pattern_speed + sigma * pattern_speed_up)
        diff_up_down = (om[i] + sigma * om_up[i]) - (pattern_speed - sigma * pattern_speed_down)
        diff_down_up = (om[i] - sigma * om_down[i]) - (pattern_speed + sigma * pattern_speed_up)
        diff_down_down = (om[i] - sigma * om_down[i]) - (pattern_speed - sigma * pattern_speed_down)

        if 0 < np.sum([diff_up_up > 0, diff_up_down > 0, diff_down_up > 0, diff_down_down > 0]) < 4:

            om_corot.append(om[i])
            r_corot.append(r[i])

            if i == 0:
                r_lower = 0
            else:
                r_lower = r[i] - (r[i] - r[i - 1]) / 2

            if i == len(om) - 1:
                r_upper = 2 * r[i] - r_lower
            else:
                r_upper = (r[i] + r[i + 1]) / 2

            if len(corotation.keys()) == 0:
                corotation[idx] = [[i], [r_lower, r_upper]]
                idx += 1
                continue

            latest_key = list(corotation.keys())[-1]

            if (i - 1) in corotation[latest_key][0]:

                corotation[latest_key][0].append(i)
                corotation[latest_key][1][-1] = r_upper

            else:

                corotation[idx] = [[i], [r_lower, r_upper]]
                idx += 1

    return om_corot, r_corot, corotation


matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

emission_to_use = ['alma_co', 'muse_mass', 'muse_ha']

# Assume the various parameters are the first in the list

star_mask = star_masks[0]
mask_outside_bar = mask_outside_bars[0]

start = time.time()

sigma = 2  # RMS threshold for associating and classifying corotation radii

rot_table = Table.read('rotation_curves/RCtable_Nov2019.fits')

environment_table = np.loadtxt('environment/PHANGSmasks_v2.dat',
                               skiprows=4,
                               dtype=str,
                               usecols=(0, 10, 11, 12, 13))

# galaxies = ['NGC1365']
# galaxies = ['IC1954']

for galaxy in alma_galaxies:

    galaxy = galaxy.strip()

    print('Plotting %s' % galaxy)

    has_alma = False
    has_muse_mass = False
    has_muse_ha = False
    has_alma_corotation = False
    has_muse_mass_corotation = False
    has_muse_ha_corotation = False
    has_bar = False

    if galaxy == 'NGC0628':
        galaxy_edit = 'NGC628'
    else:
        galaxy_edit = galaxy

    # Pull the rotation curve out of the table

    rot_idx = np.where(rot_table['Galaxy'] == galaxy.lower())

    if len(rot_idx[0]) == 0:
        print('Rotation curve not fitted')
        continue

    r = rot_table['Radius'][rot_idx]
    v = rot_table['Vrot'][rot_idx]
    v_up = rot_table['Vrot_upper'][rot_idx]
    v_down = rot_table['Vrot_lower'][rot_idx]

    om = v / r
    om_up = v_up / r
    om_down = v_down / r

    # Pull in the distance and inclination

    for i in range(len(galaxy_table)):
        if galaxy_table['NAME'][i].strip() == galaxy:
            master_idx = i
            break

    dist = galaxy_table['DIST'][master_idx]
    inc = galaxy_table['ORIENT_INCL'][master_idx]
    pa = galaxy_table['ORIENT_POSANG'][master_idx]
    r_25 = galaxy_table['SIZE_OPT_R25'][master_idx]

    r_25 *= dist * np.sin(1 / 3600 * np.pi / 180) * 1000

    # Pull in the bar radius

    try:

        environ_idx = np.where(environment_table[:, 0] == galaxy)
        bar_pa, bar_r, bar_flag = float(environment_table[environ_idx, 2][0][0]), \
                                  float(environment_table[environ_idx, 3][0][0]), \
                                  int(float(environment_table[environ_idx, 4][0][0]))

        # Take only places where bars exist, and only the Most Quality bars

        # TODO: Take quality 2 bars too (but plot different colour)

        if bar_pa > -999 and bar_flag == 1:
            # Convert R to physical scale
            bar_r *= dist * np.sin(1 / 3600 * np.pi / 180) * 1000

            # Deproject the bar

            x_rot = bar_r * np.cos(np.radians(bar_pa - pa)) * np.sin(np.radians(bar_pa - pa)) + \
                    bar_r * np.sin(np.radians(bar_pa - pa)) * np.cos(np.radians(bar_pa - pa))

            x_rot /= np.cos(np.radians(inc))
            y_rot = -bar_r * np.cos(np.radians(bar_pa - pa)) * np.cos(np.radians(bar_pa - pa)) + \
                    bar_r * np.sin(np.radians(bar_pa - pa)) * np.sin(np.radians(bar_pa - pa))

            bar_r = np.sqrt(x_rot ** 2 + y_rot ** 2)
            has_bar = True

    except IndexError:

        pass

    # Pull in the pattern speeds for MUSE mass/Halpha and ALMA (if they exist)

    try:

        file_name = muse_output + muse_version + '/mass_'

        if star_mask:
            file_name += 'smask_'
        else:
            file_name += 'nosmask_'

        if mask_outside_bar:
            file_name += 'bmask/'
        else:
            file_name += 'nobmask/'

        file_name += galaxy_edit + '_mass_'

        if star_mask:
            file_name += 'smask_'
        else:
            file_name += 'nosmask_'

        if mask_outside_bar:
            file_name += 'bmask_'
        else:
            file_name += 'nobmask_'

        file_name += 'pattern_speed_muse.txt'

        om_muse_mass, om_muse_mass_up, om_muse_mass_down = np.loadtxt(file_name, unpack=True)
        has_muse_mass = True

        # Find the corotation radii (and errors) associated with this

        has_muse_mass_corotation = False

        om_corot_muse_mass, r_corot_muse_mass, corotation_muse_mass = find_corotation_radii(om, om_up, om_down,
                                                                                            om_muse_mass,
                                                                                            om_muse_mass_up,
                                                                                            om_muse_mass_down,
                                                                                            sigma=sigma)

        if len(corotation_muse_mass.keys()) > 0:
            has_muse_mass_corotation = True

            # Calculate the corotations and errors, and write out

            corotation_muse_mass_means = []
            corotation_muse_mass_errs = []

            for key in corotation_muse_mass.keys():
                corotation_muse_mass_means.append(np.mean(corotation_muse_mass[key][-1]))
                corotation_muse_mass_errs.append(corotation_muse_mass[key][-1][-1] - corotation_muse_mass_means[-1])

            np.savetxt(corot_output + galaxy + '_muse_mass_corot.txt',
                       np.c_[corotation_muse_mass_means, corotation_muse_mass_errs],
                       header='r_corot, err_r_corot (all km/s/kpc)')

    except (FileNotFoundError, OSError):
        print('No MUSE mass found')
        pass

    # MUSE Halpha

    try:

        file_name = muse_output + muse_version + '/ha_'

        if star_mask:
            file_name += 'smask_'
        else:
            file_name += 'nosmask_'

        if mask_outside_bar:
            file_name += 'bmask/'
        else:
            file_name += 'nobmask/'

        file_name += galaxy_edit + '_ha_'

        if star_mask:
            file_name += 'smask_'
        else:
            file_name += 'nosmask_'

        if mask_outside_bar:
            file_name += 'bmask_'
        else:
            file_name += 'nobmask_'

        file_name += 'pattern_speed_muse.txt'

        om_muse_ha, om_muse_ha_up, om_muse_ha_down = np.loadtxt(file_name, unpack=True)
        has_muse_ha = True

        # Find the corotation radii (and errors) associated with this

        has_muse_ha_corotation = False

        om_corot_muse_ha, r_corot_muse_ha, corotation_muse_ha = find_corotation_radii(om, om_up, om_down,
                                                                                      om_muse_ha,
                                                                                      om_muse_ha_up,
                                                                                      om_muse_ha_down,
                                                                                      sigma=sigma)

        if len(corotation_muse_ha.keys()) > 0:
            has_muse_ha_corotation = True

            # Calculate the corotations and errors, and write out

            corotation_muse_ha_means = []
            corotation_muse_ha_errs = []

            for key in corotation_muse_ha.keys():
                corotation_muse_ha_means.append(np.mean(corotation_muse_ha[key][-1]))
                corotation_muse_ha_errs.append(corotation_muse_ha[key][-1][-1] - corotation_muse_ha_means[-1])

            np.savetxt(corot_output + galaxy + '_muse_ha_corot.txt',
                       np.c_[corotation_muse_ha_means, corotation_muse_ha_errs],
                       header='r_corot, err_r_corot (all km/s/kpc)')

    except (FileNotFoundError, OSError):
        print('No MUSE Ha found')
        pass

    try:

        file_name = alma_output + alma_version + '/'

        if mask_outside_bar:
            file_name += 'bmask/'
        else:
            file_name += 'nobmask/'

        file_name += galaxy + '_'

        if mask_outside_bar:
            file_name += 'bmask_'
        else:
            file_name += 'nobmask_'

        file_name += 'pattern_speed_alma.txt'

        om_alma, om_alma_up, om_alma_down = np.loadtxt(file_name, unpack=True)
        has_alma = True

        # Find the corotation radii (and errors) associated with this

        om_corot_alma, r_corot_alma, corotation_alma = find_corotation_radii(om, om_up, om_down,
                                                                             om_alma, om_alma_up, om_alma_down,
                                                                             sigma=sigma)

        if len(corotation_alma.keys()) > 0:
            has_alma_corotation = True

            # Calculate the corotations and errors, and write out

            corotation_alma_means = []
            corotation_alma_errs = []

            for key in corotation_alma.keys():
                corotation_alma_means.append(np.mean(corotation_alma[key][-1]))
                corotation_alma_errs.append(corotation_alma[key][-1][-1] - corotation_alma_means[-1])

            np.savetxt(corot_output + galaxy + '_alma_corot.txt',
                       np.c_[corotation_alma_means, corotation_alma_errs],
                       header='r_corot, err_r_corot (all km/s/kpc)')

    except (FileNotFoundError, OSError):
        print('No ALMA found')
        pass

    if not has_alma and not has_muse_mass and not has_muse_ha:
        print('No pattern speed found')

    fig, ax = plt.subplots(figsize=(8, 6))

    if has_muse_mass:
        plt.axhline(om_muse_mass, c='r', label=r'MUSE $M_\ast$')
        plt.axhline(om_muse_mass + sigma * om_muse_mass_up, c='r', ls='--')
        plt.axhline(om_muse_mass - sigma * om_muse_mass_down, c='r', ls='--')

    if has_muse_ha:
        plt.axhline(om_muse_ha, c='cyan', label=r'MUSE H$\alpha$')
        plt.axhline(om_muse_ha + sigma * om_muse_ha_up, c='cyan', ls='--')
        plt.axhline(om_muse_ha - sigma * om_muse_ha_down, c='cyan', ls='--')

    if has_alma:
        plt.axhline(om_alma, c='b', label='ALMA')
        plt.axhline(om_alma + sigma * om_alma_up, c='b', ls='--')
        plt.axhline(om_alma - sigma * om_alma_down, c='b', ls='--')

    if has_bar:
        plt.axvline(bar_r, c='k', ls='--', label=r'$R_\mathrm{bar}$')

    plt.errorbar(r, om, yerr=[sigma * om_down, sigma * om_up], ls='none', marker='o', c='k')

    if has_alma_corotation:
        plt.scatter(r_corot_alma, om_corot_alma, c='b', zorder=99)

        for key in corotation_alma.keys():
            plt.axvline(np.mean(corotation_alma[key][-1]), c='b')

            for value in corotation_alma[key][-1]:
                plt.axvline(value, c='b', ls='--')

    if has_muse_mass_corotation:
        plt.scatter(r_corot_muse_mass, om_corot_muse_mass, c='r', zorder=99)

        for key in corotation_muse_mass.keys():

            plt.axvline(np.mean(corotation_muse_mass[key][-1]), c='r')

            for value in corotation_muse_mass[key][-1]:
                plt.axvline(value, c='r', ls='--')

    if has_muse_ha_corotation:
        plt.scatter(r_corot_muse_ha, om_corot_muse_ha, c='cyan', zorder=99)

        for key in corotation_muse_ha.keys():

            plt.axvline(np.mean(corotation_muse_ha[key][-1]), c='cyan')

            for value in corotation_muse_ha[key][-1]:
                plt.axvline(value, c='cyan', ls='--')

    plt.ylim([0, 200])

    plt.xlabel(r'$R$ (kpc)')
    plt.ylabel(r'$\Omega$ (km s$^{-1}$ kpc$^{-1}$)')
    plt.legend(loc='upper right', frameon=False)

    ax_2 = ax.secondary_xaxis('top', functions=(lambda r: r / r_25, lambda r: r / r_25))
    ax_2.set_xlabel(r'$R/R_{25, \mathrm{opt}}$')

    # plt.show()

    plt.savefig(corot_plot + galaxy + '_corot.png',
                bbox_inches='tight')
    plt.savefig(corot_plot + galaxy + '_corot.pdf',
                bbox_inches='tight')

    plt.close()

print('Complete! Took %.2fs' % (time.time() - start))
