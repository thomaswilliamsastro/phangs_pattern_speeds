# -*- coding: utf-8 -*-
"""
Pull everything together into a .fits table (and LaTeX table)

@author: Tom Williams
"""

import os

import astropy.units as u
import numpy as np
from astropy.io import ascii
from astropy.table import Table

from vars import phangs_folder, alma_version, muse_version, corot_version, alma_galaxies, output_folder, \
    alma_output, muse_output, corot_output, galaxy_table, star_masks, mask_outside_bars

emission_to_include = ['alma_co', 'muse_mass', 'muse_ha']

# Assume the various parameters are the first in the list

star_mask = star_masks[0]
mask_outside_bar = mask_outside_bars[0]

print('Creating PHANGS corotation table %s' % corot_version)
print('MUSE %s:\n-star mask: %s\n-beyond bar mask: %s'
      % (muse_version, star_mask, mask_outside_bar))
print('ALMA %s:\n-beyond bar mask: %s'
      % (alma_version, mask_outside_bar))

os.chdir(phangs_folder)

n_alma_corotations = 0
n_muse_corotations = 0

keys = None
alma_corot_keys = []
muse_mass_corot_keys = []
muse_ha_corot_keys = []
data = {}

# alma_galaxies = ['NGC0628   ']

alma_co_q_table = ascii.read(output_folder + 'qf_alma_co.csv')
muse_mass_q_table = ascii.read(output_folder + 'qf_muse_mass.csv')
muse_ha_q_table = ascii.read(output_folder + 'qf_muse_ha.csv')

for galaxy in alma_galaxies:
    has_alma_pattern_speed = False
    has_muse_mass_pattern_speed = False
    has_muse_ha_pattern_speed = False

    galaxy_row = galaxy_table[galaxy_table['NAME'] == galaxy]

    # Pull out distance and inclination from the master table, so users can convert back if necessary
    dist = galaxy_row['DIST'][0]
    incl = galaxy_row['ORIENT_INCL'][0]
    pgc = galaxy_row['PGC'][0]

    galaxy = galaxy.strip()
    # print(galaxy)

    # Pull out the quality flags.
    alma_co_row = alma_co_q_table[alma_co_q_table['GALAXY'] == galaxy]
    muse_mass_row = muse_mass_q_table[muse_mass_q_table['GALAXY'] == galaxy]
    muse_ha_row = muse_ha_q_table[muse_ha_q_table['GALAXY'] == galaxy]

    # TODO: At some point this will be the mode, for now just use my flag.
    try:
        alma_co_q = alma_co_row['FLAG_TW'][0]
    except IndexError:
        alma_co_q = np.nan
    try:
        muse_mass_q = muse_mass_row['FLAG_TW'][0]
    except IndexError:
        muse_mass_q = np.nan
    try:
        muse_ha_q = muse_ha_row['FLAG_TW'][0]
    except IndexError:
        muse_ha_q = np.nan

    # Read in the ALMA pattern speed

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
        has_alma_pattern_speed = True
    except OSError:
        om_alma, om_alma_up, om_alma_down = np.nan, np.nan, np.nan

    # Read in the MUSE stellar mass pattern speed

    try:

        if galaxy == 'NGC0628':
            galaxy_edit = 'NGC628'
        else:
            galaxy_edit = galaxy

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

        has_muse_mass_pattern_speed = True
    except OSError:
        om_muse_mass, om_muse_mass_up, om_muse_mass_down = np.nan, np.nan, np.nan

    # Read in the MUSE Halpha pattern speed

    try:

        if galaxy == 'NGC0628':
            galaxy_edit = 'NGC628'
        else:
            galaxy_edit = galaxy

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

        has_muse_ha_pattern_speed = True
    except OSError:
        om_muse_ha, om_muse_ha_up, om_muse_ha_down = np.nan, np.nan, np.nan

    if not has_alma_pattern_speed and not has_muse_mass_pattern_speed and not has_muse_ha_pattern_speed:
        continue

    try:
        data['GALAXY'].append(galaxy)
    except KeyError:
        data['GALAXY'] = [galaxy]
    try:
        data['PGC'].append(pgc)
    except KeyError:
        data['PGC'] = [pgc]
    try:
        data['DIST'].append(dist)
    except KeyError:
        data['DIST'] = [dist]
    try:
        data['INCL'].append(incl)
    except KeyError:
        data['INCL'] = [incl]

    try:
        data['OM_P_ALMA'].append(om_alma)
    except KeyError:
        data['OM_P_ALMA'] = [om_alma]
    try:
        data['OM_P_ALMA_ERR_UP'].append(om_alma_up)
    except KeyError:
        data['OM_P_ALMA_ERR_UP'] = [om_alma_up]
    try:
        data['OM_P_ALMA_ERR_DOWN'].append(om_alma_down)
    except KeyError:
        data['OM_P_ALMA_ERR_DOWN'] = [om_alma_down]
    try:
        data['OM_P_ALMA_QUAL'].append(alma_co_q)
    except KeyError:
        data['OM_P_ALMA_QUAL'] = [alma_co_q]

    try:
        data['OM_P_MUSE_MASS'].append(om_muse_mass)
    except KeyError:
        data['OM_P_MUSE_MASS'] = [om_muse_mass]
    try:
        data['OM_P_MUSE_MASS_ERR_UP'].append(om_muse_mass_up)
    except KeyError:
        data['OM_P_MUSE_MASS_ERR_UP'] = [om_muse_mass_up]
    try:
        data['OM_P_MUSE_MASS_ERR_DOWN'].append(om_muse_mass_down)
    except KeyError:
        data['OM_P_MUSE_MASS_ERR_DOWN'] = [om_muse_mass_down]
    try:
        data['OM_P_MUSE_MASS_QUAL'].append(muse_mass_q)
    except KeyError:
        data['OM_P_MUSE_MASS_QUAL'] = [muse_mass_q]

    try:
        data['OM_P_MUSE_HA'].append(om_muse_ha)
    except KeyError:
        data['OM_P_MUSE_HA'] = [om_muse_ha]
    try:
        data['OM_P_MUSE_HA_ERR_UP'].append(om_muse_ha_up)
    except KeyError:
        data['OM_P_MUSE_HA_ERR_UP'] = [om_muse_ha_up]
    try:
        data['OM_P_MUSE_HA_ERR_DOWN'].append(om_muse_ha_down)
    except KeyError:
        data['OM_P_MUSE_HA_ERR_DOWN'] = [om_muse_ha_down]
    try:
        data['OM_P_MUSE_HA_QUAL'].append(muse_ha_q)
    except KeyError:
        data['OM_P_MUSE_HA_QUAL'] = [muse_ha_q]

    if not keys:
        keys = list(data.keys())

    # Check for any ALMA corotation radii

    try:
        alma_corot, alma_corot_err = np.loadtxt(corot_output + galaxy + '_alma_corot.txt',
                                                unpack=True)

        # Loop over the corotations and put into dictionary as necessary

        alma_cols_included = []

        if type(alma_corot) == np.ndarray:

            for i in range(len(alma_corot)):

                corot_colname = 'R_CR_ALMA_%d' % (i + 1)
                corot_err_colname = corot_colname + '_ERR'

                alma_cols_included.extend([corot_colname, corot_err_colname])

                if corot_colname not in alma_corot_keys:
                    alma_corot_keys.extend([corot_colname, corot_err_colname])
                    new_col = np.zeros(len(data['GALAXY']))
                    new_col[new_col == 0] = np.nan
                    new_col = list(new_col)
                    new_err_col = new_col.copy()

                    new_col[-1] = alma_corot[i]
                    new_err_col[-1] = alma_corot_err[i]

                    data[corot_colname] = new_col
                    data[corot_err_colname] = new_err_col

                else:

                    data[corot_colname].append(alma_corot[i])
                    data[corot_err_colname].append(alma_corot_err[i])

        else:

            corot_colname = 'R_CR_ALMA_1'
            corot_err_colname = corot_colname + '_ERR'

            alma_cols_included.extend([corot_colname, corot_err_colname])

            if corot_colname not in alma_corot_keys:
                alma_corot_keys.extend([corot_colname, corot_err_colname])
                new_col = np.zeros(len(data['GALAXY']))
                new_col[new_col == 0] = np.nan
                new_col = list(new_col)
                new_err_col = new_col.copy()

                new_col[-1] = alma_corot
                new_err_col[-1] = alma_corot_err

                data[corot_colname] = new_col
                data[corot_err_colname] = new_err_col

            else:

                data[corot_colname].append(alma_corot)
                data[corot_err_colname].append(alma_corot_err)

        for key in alma_corot_keys:

            if key not in alma_cols_included:
                data[key].append(np.nan)

    except OSError:

        for key in alma_corot_keys:
            data[key].append(np.nan)

    # Check for any MUSE mass corotation radii

    try:

        muse_mass_corot, muse_mass_corot_err = np.loadtxt(corot_output + galaxy + '_muse_mass_corot.txt',
                                                          unpack=True)

        # Loop over the corotations and put into dictionary as necessary

        muse_mass_cols_included = []

        if type(muse_mass_corot) == np.ndarray:

            for i in range(len(muse_mass_corot)):

                corot_colname = 'R_CR_MUSE_MASS_%d' % (i + 1)
                corot_err_colname = corot_colname + '_ERR'

                muse_mass_cols_included.extend([corot_colname, corot_err_colname])

                if corot_colname not in muse_mass_corot_keys:
                    muse_mass_corot_keys.extend([corot_colname, corot_err_colname])
                    new_col = np.zeros(len(data['GALAXY']))
                    new_col[new_col == 0] = np.nan
                    new_col = list(new_col)
                    new_err_col = new_col.copy()

                    new_col[-1] = muse_mass_corot[i]
                    new_err_col[-1] = muse_mass_corot_err[i]

                    data[corot_colname] = new_col
                    data[corot_err_colname] = new_err_col

                else:

                    data[corot_colname].append(muse_mass_corot[i])
                    data[corot_err_colname].append(muse_mass_corot_err[i])

        else:

            corot_colname = 'R_CR_MUSE_MASS_1'
            corot_err_colname = corot_colname + '_ERR'

            muse_mass_cols_included.extend([corot_colname, corot_err_colname])

            if corot_colname not in muse_mass_corot_keys:
                muse_mass_corot_keys.extend([corot_colname, corot_err_colname])
                new_col = np.zeros(len(data['GALAXY']))
                new_col[new_col == 0] = np.nan
                new_col = list(new_col)
                new_err_col = new_col.copy()

                new_col[-1] = muse_mass_corot
                new_err_col[-1] = muse_mass_corot_err

                data[corot_colname] = new_col
                data[corot_err_colname] = new_err_col

            else:

                data[corot_colname].append(muse_mass_corot)
                data[corot_err_colname].append(muse_mass_corot_err)

        for key in muse_mass_corot_keys:

            if key not in muse_mass_cols_included:
                data[key].append(np.nan)

    except OSError:

        for key in muse_mass_corot_keys:
            data[key].append(np.nan)

    # Check for any MUSE Ha corotation radii

    try:

        muse_ha_corot, muse_ha_corot_err = np.loadtxt(corot_output + galaxy + '_muse_ha_corot.txt',
                                                      unpack=True)

        # Loop over the corotations and put into dictionary as necessary

        muse_ha_cols_included = []

        if type(muse_ha_corot) == np.ndarray:

            for i in range(len(muse_ha_corot)):

                corot_colname = 'R_CR_MUSE_HA_%d' % (i + 1)
                corot_err_colname = corot_colname + '_ERR'

                muse_ha_cols_included.extend([corot_colname, corot_err_colname])

                if corot_colname not in muse_ha_corot_keys:
                    muse_ha_corot_keys.extend([corot_colname, corot_err_colname])
                    new_col = np.zeros(len(data['GALAXY']))
                    new_col[new_col == 0] = np.nan
                    new_col = list(new_col)
                    new_err_col = new_col.copy()

                    new_col[-1] = muse_ha_corot[i]
                    new_err_col[-1] = muse_ha_corot_err[i]

                    data[corot_colname] = new_col
                    data[corot_err_colname] = new_err_col

                else:

                    data[corot_colname].append(muse_ha_corot[i])
                    data[corot_err_colname].append(muse_ha_corot_err[i])

        else:

            corot_colname = 'R_CR_MUSE_HA_1'
            corot_err_colname = corot_colname + '_ERR'

            muse_ha_cols_included.extend([corot_colname, corot_err_colname])

            if corot_colname not in muse_ha_corot_keys:
                muse_ha_corot_keys.extend([corot_colname, corot_err_colname])
                new_col = np.zeros(len(data['GALAXY']))
                new_col[new_col == 0] = np.nan
                new_col = list(new_col)
                new_err_col = new_col.copy()

                new_col[-1] = muse_ha_corot
                new_err_col[-1] = muse_ha_corot_err

                data[corot_colname] = new_col
                data[corot_err_colname] = new_err_col

            else:

                data[corot_colname].append(muse_ha_corot)
                data[corot_err_colname].append(muse_ha_corot_err)

        for key in muse_ha_corot_keys:

            if key not in muse_ha_cols_included:
                data[key].append(np.nan)

    except OSError:

        for key in muse_ha_corot_keys:
            data[key].append(np.nan)

keys.extend(alma_corot_keys)
keys.extend(muse_mass_corot_keys)
keys.extend(muse_ha_corot_keys)

table = Table()

for key in keys:

    column = data[key]

    if 'DIST' in key:
        column = column * u.Mpc
    elif 'INCL' in key:
        column = column * u.degree
    elif 'OM_P' in key and 'QUAL' not in key:
        column = column * u.km / u.s / u.kpc
    elif 'R_CR' in key:
        column = column * u.kpc

    table.add_column(column, name=key)

table.write(output_folder + 'pattern_speed_table_' + corot_version + '.fits',
            format='fits', overwrite=True)

# Secondly, TeXiFy this table for public consumption

names = ['Galaxy', 'PGC', 'D', '$i$']
units = {'D': u.Mpc,
         '$i$': u.deg,
         r'$\Omega_{\rm P, ALMA}$': u.km * u.s ** -1 * u.kpc ** -1,
         r'$\Omega_{\rm P, MUSE mass}$': u.km * u.s ** -1 * u.kpc ** -1,
         r'$\Omega_{\rm P, MUSE H-alpha}$': u.km * u.s ** -1 * u.kpc ** -1,
         }

latex_table = Table()
latex_table.add_column(table['GALAXY'], name='Galaxy')
latex_table.add_column(table['PGC'], name='PGC')
latex_table.add_column(table['DIST'], name='D')
latex_table.add_column(table['INCL'], name='i')

latex_alma_om_p = ['$%.1f^{+%.1f}_{-%.1f}$' % (table['OM_P_ALMA'][i],
                                               table['OM_P_ALMA_ERR_UP'][i],
                                               table['OM_P_ALMA_ERR_DOWN'][i])
                   for i in range(len(table))]

for i in range(len(latex_alma_om_p)):
    if 'nan' in latex_alma_om_p[i]:
        latex_alma_om_p[i] = 'nan'

latex_table.add_column(latex_alma_om_p, name='om_p_alma')
latex_table.add_column(table['OM_P_ALMA_QUAL'].astype(int), name='om_p_alma_qual')

names.extend([r'$\Omega_{\rm P, ALMA}$', r'Q$_{\rm ALMA}$'])

# Include the ALMA corotations

alma_corot_names = []
alma_corot_number = 1

for key in alma_corot_keys:

    if '_ERR' not in key:
        new_col = [r'$%.1f\pm%.1f$' % (table[key][i],
                                       table[key + '_ERR'][i])
                   for i in range(len(table))]

        for i in range(len(new_col)):
            if 'nan' in new_col[i]:
                new_col[i] = 'nan'

        latex_table.add_column(new_col, name=key.lower())
        new_name = r'$R_{\rm CR%d, A}$' % alma_corot_number
        alma_corot_names.append(new_name)

        units.update({new_name: u.kpc})

        alma_corot_number += 1

names.extend(alma_corot_names)

latex_muse_mass_om_p = ['$%.1f^{+%.1f}_{-%.1f}$' % (table['OM_P_MUSE_MASS'][i],
                                                    table['OM_P_MUSE_MASS_ERR_UP'][i],
                                                    table['OM_P_MUSE_MASS_ERR_DOWN'][i])
                        for i in range(len(table))]

for i in range(len(latex_muse_mass_om_p)):
    if 'nan' in latex_muse_mass_om_p[i]:
        latex_muse_mass_om_p[i] = 'nan'
latex_table.add_column(latex_muse_mass_om_p, name='om_p_muse_mass')
latex_table.add_column(table['OM_P_MUSE_MASS_QUAL'].astype(int), name='om_p_muse_mass_qual')

names.extend([r'$\Omega_{\rm P, MUSE mass}$', r'Q$_{\rm MUSE mass}$'])

# Include the MUSE corotations

muse_mass_corot_names = []
muse_mass_corot_number = 1

for key in muse_mass_corot_keys:

    if '_ERR' not in key:
        new_col = [r'$%.1f\pm%.1f$' % (table[key][i],
                                       table[key + '_ERR'][i])
                   for i in range(len(table))]

        for i in range(len(new_col)):
            if 'nan' in new_col[i]:
                new_col[i] = 'nan'

        latex_table.add_column(new_col, name=key.lower())
        new_name = r'$R_{\rm CR%d, M M}$' % muse_mass_corot_number
        muse_mass_corot_names.append(new_name)

        units.update({new_name: u.kpc})

        muse_mass_corot_number += 1

names.extend(muse_mass_corot_names)

latex_muse_ha_om_p = ['$%.1f^{+%.1f}_{-%.1f}$' % (table['OM_P_MUSE_HA'][i],
                                                  table['OM_P_MUSE_HA_ERR_UP'][i],
                                                  table['OM_P_MUSE_HA_ERR_DOWN'][i])
                      for i in range(len(table))]

for i in range(len(latex_muse_ha_om_p)):
    if 'nan' in latex_muse_ha_om_p[i]:
        latex_muse_ha_om_p[i] = 'nan'
latex_table.add_column(latex_muse_ha_om_p, name='om_p_muse_ha')
latex_table.add_column(table['OM_P_MUSE_HA_QUAL'].astype(int), name='om_p_muse_ha_qual')

names.extend([r'$\Omega_{\rm P, MUSE Ha}$', r'Q$_{\rm MUSE Ha}$'])

# Include the MUSE corotations

muse_ha_corot_names = []
muse_ha_corot_number = 1

for key in muse_ha_corot_keys:

    if '_ERR' not in key:
        new_col = [r'$%.1f\pm%.1f$' % (table[key][i],
                                       table[key + '_ERR'][i])
                   for i in range(len(table))]

        for i in range(len(new_col)):
            if 'nan' in new_col[i]:
                new_col[i] = 'nan'

        latex_table.add_column(new_col, name=key.lower())
        new_name = r'$R_{\rm CR%d, M H}$' % muse_ha_corot_number
        muse_ha_corot_names.append(new_name)

        units.update({new_name: u.kpc})

        muse_ha_corot_number += 1

names.extend(muse_ha_corot_names)

for key in units.keys():
    units[key] = units[key].to_string('latex')

ascii.write(latex_table,
            output_folder + 'pattern_speed_table_' + corot_version + '.tex',
            # sys.stdout,
            names=names,
            fill_values=[('nan', r'\nodata'),
                         ('-9223372036854775808', r'\nodata')],
            format='aastex',
            overwrite=True,
            latexdict={'units': units,
                       'tabletype': 'deluxetable*',
                       'preamble': r'\tablecaption{Pattern speeds and calculated co-rotation radii for 80 galaxies}'
                                   '\n\label{table:pattern_speeds}'},
            )

print('Complete!')
