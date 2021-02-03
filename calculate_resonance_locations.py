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

from vars import phangs_folder, alma_version, muse_version, alma_output, muse_output, galaxy_table, \
    star_masks, mask_outside_bars, resonance_output, resonance_plot


def find_resonance_radii(om, om_up, om_down, r, pattern_speed, pattern_speed_up, pattern_speed_down,
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


def write_resonances(resonance_dict, file_name):
    """ Calculate the corotations and errors, and write out.

    """

    resonance_means = []
    resonance_errs = []

    for key in resonance_dict.keys():
        resonance_means.append(np.mean(resonance_dict[key][-1]))
        resonance_errs.append(resonance_dict[key][-1][-1] - resonance_means[-1])

    np.savetxt(file_name,
               np.c_[resonance_means, resonance_errs],
               header='r_resonance, err_r_resonance (all km/s/kpc)')

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
old_sample_table = Table.read('documents/phangs_sample_table_v1p4.fits')

if not os.path.exists(resonance_output):
    os.makedirs(resonance_output)
if not os.path.exists(resonance_plot):
    os.makedirs(resonance_plot)

alma_galaxies = ['ngc3351']

for galaxy in alma_galaxies:

    galaxy = galaxy.upper()

    print('Plotting %s' % galaxy)

    has_alma = False
    has_muse_mass = False
    has_muse_ha = False
    has_alma_corotation = False
    has_muse_mass_corotation = False
    has_muse_ha_corotation = False
    has_bar = False

    # Pull in the distance and inclination

    galaxy_row = galaxy_table[galaxy_table['name'] == galaxy.lower()]

    dist = galaxy_row['dist'][0]
    inc = galaxy_row['orient_incl'][0]
    pa = galaxy_row['orient_posang'][0]
    r_25 = galaxy_row['size_r25'][0]

    r_25 *= dist * np.sin(1 / 3600 * np.pi / 180) * 1000

    # Pull the rotation curve out of the table

    rot_idx = np.where(rot_table['Galaxy'] == galaxy.lower())

    if len(rot_idx[0]) == 0:
        print('Rotation curve not fitted')
        continue

    r = rot_table['Radius'][rot_idx]
    v = rot_table['Vrot'][rot_idx]
    v_up = rot_table['Vrot_upper'][rot_idx]
    v_down = rot_table['Vrot_lower'][rot_idx]

    # We need to edit the R to match the newest distances

    old_sample_table_idx = np.where(old_sample_table['NAME'] == galaxy.upper().ljust(10, ' '))[0][0]
    old_dist = old_sample_table['DIST'][old_sample_table_idx]

    r *= dist/old_dist

    om = v / r
    om_up = v_up / r
    om_down = v_down / r

    # And splines

    spline_model = Table.read('rotation_curve_spline/' + galaxy + '_spline_model.ecsv', format='ascii.ecsv')
    v_spline = spline_model['v_circ']
    r_spline = spline_model['r_gal'] * dist * np.sin(1 / 3600 * np.pi / 180) * 1000

    dr = r_spline[1] - r_spline[0]
    om_spline = v_spline / r_spline
    kappa_spline = np.sqrt((2 * om_spline / r_spline) * np.gradient(r_spline ** 2 * om_spline, dr))

    ilr_spline = om_spline - kappa_spline / 2
    olr_spline = om_spline + kappa_spline / 2

    lr_err = np.zeros_like(ilr_spline)

    # Pull in the bar radius

    bar_pa, bar_r, bar_flag = galaxy_row['morph_bar_pa'][0], galaxy_row['morph_bar_r'][0], \
                              galaxy_row['morph_bar_flag'][0]

    # Take only places where bars exist, and only the Most Quality bars

    if ~np.isnan(bar_r) and bar_flag == 1:
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

        file_name += galaxy + '_mass_'

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

        om_corot_muse_mass, r_corot_muse_mass, corotation_muse_mass = find_resonance_radii(om, om_up, om_down, r,
                                                                                           om_muse_mass,
                                                                                           om_muse_mass_up,
                                                                                           om_muse_mass_down,
                                                                                           sigma=sigma)

        if len(corotation_muse_mass.keys()) > 0:
            has_muse_mass_corotation = True

            write_resonances(corotation_muse_mass, resonance_output + galaxy + '_muse_mass_corot.txt')

        # Do the same for the ILR/OLR

        om_ilr_muse_mass, r_ilr_muse_mass, ilr_muse_mass = find_resonance_radii(ilr_spline, lr_err, lr_err,
                                                                                r_spline,
                                                                                om_muse_mass,
                                                                                om_muse_mass_up,
                                                                                om_muse_mass_down,
                                                                                sigma=sigma)

        if len(ilr_muse_mass.keys()) > 0:
            write_resonances(ilr_muse_mass, resonance_output + galaxy + '_muse_mass_ilr.txt')

        om_olr_muse_mass, r_olr_muse_mass, olr_muse_mass = find_resonance_radii(olr_spline, lr_err, lr_err,
                                                                                r_spline,
                                                                                om_muse_mass,
                                                                                om_muse_mass_up,
                                                                                om_muse_mass_down,
                                                                                sigma=sigma)

        if len(olr_muse_mass.keys()) > 0:
            write_resonances(olr_muse_mass, resonance_output + galaxy + '_muse_mass_olr.txt')

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

        file_name += galaxy + '_ha_'

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

        om_corot_muse_ha, r_corot_muse_ha, corotation_muse_ha = find_resonance_radii(om, om_up, om_down, r,
                                                                                     om_muse_ha,
                                                                                     om_muse_ha_up,
                                                                                     om_muse_ha_down,
                                                                                     sigma=sigma)

        if len(corotation_muse_ha.keys()) > 0:
            has_muse_ha_corotation = True

            write_resonances(corotation_muse_ha, resonance_output + galaxy + '_muse_ha_corot.txt')

        # Do the same for the ILR/OLR

        om_ilr_muse_ha, r_ilr_muse_ha, ilr_muse_ha = find_resonance_radii(ilr_spline, lr_err, lr_err,
                                                                          r_spline,
                                                                          om_muse_ha,
                                                                          om_muse_ha_up,
                                                                          om_muse_ha_down,
                                                                          sigma=sigma)

        if len(ilr_muse_ha.keys()) > 0:
            write_resonances(ilr_muse_ha, resonance_output + galaxy + '_muse_ha_ilr.txt')

        om_olr_muse_ha, r_olr_muse_ha, olr_muse_ha = find_resonance_radii(olr_spline, lr_err, lr_err,
                                                                          r_spline,
                                                                          om_muse_ha,
                                                                          om_muse_ha_up,
                                                                          om_muse_ha_down,
                                                                          sigma=sigma)

        if len(olr_muse_ha.keys()) > 0:
            write_resonances(olr_muse_ha, resonance_output + galaxy + '_muse_ha_olr.txt')

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

        om_corot_alma, r_corot_alma, corotation_alma = find_resonance_radii(om, om_up, om_down, r,
                                                                            om_alma, om_alma_up, om_alma_down,
                                                                            sigma=sigma)

        if len(corotation_alma.keys()) > 0:
            has_alma_corotation = True

            write_resonances(corotation_alma, resonance_output + galaxy + '_alma_corot.txt')

        # Do the same for the ILR/OLR

        om_ilr_alma, r_ilr_alma, ilr_alma = find_resonance_radii(ilr_spline, lr_err, lr_err,
                                                                 r_spline,
                                                                 om_alma,
                                                                 om_alma_up,
                                                                 om_alma_down,
                                                                 sigma=sigma)

        if len(ilr_alma.keys()) > 0:
            write_resonances(ilr_alma, resonance_output + galaxy + '_alma_ilr.txt')

        om_olr_alma, r_olr_alma, olr_alma = find_resonance_radii(olr_spline, lr_err, lr_err,
                                                                 r_spline,
                                                                 om_alma,
                                                                 om_alma_up,
                                                                 om_alma_down,
                                                                 sigma=sigma)

        if len(olr_alma.keys()) > 0:
            write_resonances(olr_alma, resonance_output + galaxy + '_alma_olr.txt')

    except (FileNotFoundError, OSError):
        print('No ALMA found')
        pass

    if not has_alma and not has_muse_mass and not has_muse_ha:
        print('No pattern speed found')

    fig, ax = plt.subplots(figsize=(4, 3.5))

    xmin, xmax = 0, 1.01 * np.nanmax(r)

    plt.errorbar(r, om, yerr=[sigma * om_down, sigma * om_up], ls='none', marker='o', c='k', label='Lang+ (2020)')

    # Also plot on the lines for LRs

    plt.plot(r_spline, olr_spline, c='k', ls='--', label=r'$\Omega \pm \kappa/2$')
    plt.plot(r_spline, ilr_spline, c='k', ls='--')

    if has_bar:
        plt.axvline(bar_r, c='k', ls='-.', label=r'$R_\mathrm{bar}$')

    if has_muse_mass:
        plt.axhline(om_muse_mass, c='g', label=r'$\Omega_\mathrm{P}$')
        plt.fill_between([xmin, xmax], om_muse_mass - sigma * om_muse_mass_down, om_muse_mass + sigma * om_muse_mass_up,
                         color='g', alpha=0.25)
        # plt.axhline(om_muse_mass + sigma * om_muse_mass_up, c='r', ls='--')
        # plt.axhline(om_muse_mass - sigma * om_muse_mass_down, c='r', ls='--')

    # if has_muse_ha:
    #     plt.axhline(om_muse_ha, c='cyan', label=r'MUSE H$\alpha$')
    #     plt.axhline(om_muse_ha + sigma * om_muse_ha_up, c='cyan', ls='--')
    #     plt.axhline(om_muse_ha - sigma * om_muse_ha_down, c='cyan', ls='--')

    # if has_alma:
    #     plt.axhline(om_alma, c='b', label='ALMA')
    #     plt.axhline(om_alma + sigma * om_alma_up, c='b', ls='--')
    #     plt.axhline(om_alma - sigma * om_alma_down, c='b', ls='--')

    # if has_alma_corotation:
    #     plt.scatter(r_corot_alma, om_corot_alma, c='b', zorder=99)
    #
    #     for key in corotation_alma.keys():
    #         plt.axvline(np.mean(corotation_alma[key][-1]), c='b')
    #
    #         for value in corotation_alma[key][-1]:
    #             plt.axvline(value, c='b', ls='--')

    if has_muse_mass_corotation:
        plt.scatter(r_corot_muse_mass, om_corot_muse_mass, c='r', zorder=99)

        for key in corotation_muse_mass.keys():

            plt.axvline(np.mean(corotation_muse_mass[key][-1]), c='r', label=r'$R_\mathrm{CR}$')
            plt.fill_between(corotation_muse_mass[key][-1], 0, 1000, color='r', alpha=0.3)

    # if has_muse_ha_corotation:
    #     plt.scatter(r_corot_muse_ha, om_corot_muse_ha, c='cyan', zorder=99)
    #
    #     for key in corotation_muse_ha.keys():
    #
    #         plt.axvline(np.mean(corotation_muse_ha[key][-1]), c='cyan')
    #
    #         for value in corotation_muse_ha[key][-1]:
    #             plt.axvline(value, c='cyan', ls='--')

    plt.ylim([0, 200])
    plt.xlim(xmin, xmax)

    plt.xlabel(r'$R$ (kpc)')
    plt.ylabel(r'$\Omega$ (km s$^{-1}$ kpc$^{-1}$)')
    leg = plt.legend(loc='center left', framealpha=0.75, bbox_to_anchor=(1.1, .5))
    leg.get_frame().set_linewidth(0.0)

    ax_2 = ax.secondary_xaxis('top', functions=(lambda r: r / r_25, lambda r: r / r_25))
    ax_2.set_xlabel(r'$R/R_{25, \mathrm{opt}}$')

    # plt.tight_layout()

    # plt.show()

    plt.savefig(resonance_plot + galaxy + '_resonances.png',
                bbox_inches='tight')
    plt.savefig(resonance_plot + galaxy + '_resonances.pdf',
                bbox_inches='tight')

    no

    plt.close()

print('Complete! Took %.2fs' % (time.time() - start))
