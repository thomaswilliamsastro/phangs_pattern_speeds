# -*- coding: utf-8 -*-
"""
Pull everything together into a .fits table (and LaTeX table)

@author: Tom Williams
"""

import os

import astropy.units as u
import numpy as np
from astropy.io import ascii
from astropy.table import Table, vstack
from scipy.stats import mode

from vars import phangs_folder, alma_version, muse_version, pattern_speed_version, alma_galaxies, output_folder, \
    alma_output, muse_output, resonance_output, galaxy_table, star_masks, mask_outside_bars


def add_radii(file_name, data_dict, radii_keys, radii_name='R_CR_ALMA'):
    if os.path.exists(file_name):
        radii, radii_err = np.loadtxt(file_name, unpack=True)

        # Loop over the corotations and put into dictionary as necessary

        cols_included = []

        if type(radii) == np.ndarray:

            for i in range(len(radii)):

                radii_colname = '%s_%d' % (radii_name, i + 1)
                radii_err_colname = radii_colname + '_ERR'

                cols_included.extend([radii_colname, radii_err_colname])

                if radii_colname not in radii_keys:
                    radii_keys.extend([radii_colname, radii_err_colname])
                    new_col = np.zeros(len(data_dict['GALAXY']))
                    new_col[new_col == 0] = np.nan
                    new_col = list(new_col)
                    new_err_col = new_col.copy()

                    new_col[-1] = radii[i]
                    new_err_col[-1] = radii_err[i]

                    data_dict[radii_colname] = new_col
                    data_dict[radii_err_colname] = new_err_col

                else:

                    data_dict[radii_colname].append(radii[i])
                    data[radii_err_colname].append(radii_err[i])

        else:

            radii_colname = '%s_%d' % (radii_name, 1)
            radii_err_colname = radii_colname + '_ERR'

            cols_included.extend([radii_colname, radii_err_colname])

            if radii_colname not in radii_keys:
                radii_keys.extend([radii_colname, radii_err_colname])
                new_col = np.zeros(len(data_dict['GALAXY']))
                new_col[new_col == 0] = np.nan
                new_col = list(new_col)
                new_err_col = new_col.copy()

                new_col[-1] = radii
                new_err_col[-1] = radii_err

                data_dict[radii_colname] = new_col
                data_dict[radii_err_colname] = new_err_col

            else:

                data_dict[radii_colname].append(radii)
                data_dict[radii_err_colname].append(radii_err)

        for key in radii_keys:

            if key not in cols_included:
                data_dict[key].append(np.nan)

    else:
        for key in radii_keys:
            data_dict[key].append(np.nan)

    return data_dict, radii_keys


def texify_om_p(latex_table, names, table, data_name='ALMA', short_name='A'):
    latex_om_p = ['$%.1f^{+%.1f}_{-%.1f}$' % (table['OM_P_' + data_name][i],
                                              table['OM_P_' + data_name + '_ERR_UP'][i],
                                              table['OM_P_' + data_name + '_ERR_DOWN'][i])
                  for i in range(len(latex_table))]

    for i in range(len(latex_om_p)):
        if 'nan' in latex_om_p[i]:
            latex_om_p[i] = 'nan'

    if data_name == 'MUSE_MASS':
        omega_name = 'P'
    else:
        omega_name = 'clump'

    latex_table.add_column(latex_om_p,
                           name='om_p_' + data_name.lower())
    latex_table.add_column(table['OM_P_' + data_name + '_QUAL'].astype(int),
                           name='om_p_' + data_name.lower() + '_qual')

    names.extend([r'$\Omega_{\rm '+omega_name+', ' + short_name + '}$',
                  r'Q$_{\rm ' + short_name + '}$'])

    return latex_table, names


def texify_resonances(resonance_keys, resonance_table, latex_table, col_align,
                      names, name_style, n_cols=1):
    """TeX up resonances for the table.

    """

    resonance_names = []
    resonance_number = 1

    for key in resonance_keys:

        if '_ERR' not in key:
            new_col = [r'$%.1f\pm%.1f$' % (resonance_table[key][i],
                                           resonance_table[key + '_ERR'][i])
                       for i in range(len(resonance_table))]

            for i in range(len(new_col)):
                if 'nan' in new_col[i]:
                    new_col[i] = 'nan'

            latex_table.add_column(new_col, name=key.lower())
            new_name = name_style % resonance_number
            resonance_names.append(new_name)

            units.update({new_name: u.kpc})

            if resonance_number <= n_cols:
                col_align += 'c'
            else:
                col_align += 'h'

            resonance_number += 1

    names.extend(resonance_names)

    return latex_table, col_align, names


emission_to_include = ['alma_co', 'muse_mass', 'muse_ha']

# Assume the various parameters are the first in the list

star_mask = star_masks[0]
mask_outside_bar = mask_outside_bars[0]

print('Creating PHANGS corotation table %s' % pattern_speed_version)
print('MUSE %s:\n-star mask: %s\n-beyond bar mask: %s'
      % (muse_version, star_mask, mask_outside_bar))
print('ALMA %s:\n-beyond bar mask: %s'
      % (alma_version, mask_outside_bar))

os.chdir(phangs_folder)

n_corot_cols = 1

n_alma_corotations = 0
n_muse_corotations = 0

keys = None
alma_corot_keys = []
muse_mass_corot_keys = []
muse_ha_corot_keys = []

alma_ilr_keys = []
alma_olr_keys = []

muse_mass_ilr_keys = []
muse_mass_olr_keys = []

muse_ha_ilr_keys = []
muse_ha_olr_keys = []

data = {}

quality_flaggers = ['ee', 'es', 'tw']
quality_flags = {}

#  Create quality flags

for quality_flagger in quality_flaggers:

    tab_name = os.path.join(output_folder, 'quality_flags', 'batch1_' + quality_flagger + '.txt')
    tab1 = ascii.read(tab_name, format='no_header')

    tab_name = os.path.join(output_folder, 'quality_flags', 'batch2_' + quality_flagger + '.txt')
    tab2 = ascii.read(tab_name, format='no_header')
    tab = vstack([tab1, tab2])

    for galaxy in alma_galaxies:

        galaxy = galaxy.upper()

        if galaxy.upper() + '_alma' in tab['col1']:
            idx = np.where(tab['col1'] == galaxy + '_alma')[0][0]
            try:
                quality_flags[galaxy + '_alma'].append(tab['col2'][idx] + 1)
            except KeyError:
                quality_flags[galaxy + '_alma'] = [tab['col2'][idx] + 1]

        if galaxy.upper() + '_ha' in tab['col1']:
            idx = np.where(tab['col1'] == galaxy + '_ha')[0][0]
            try:
                quality_flags[galaxy + '_muse_ha'].append(tab['col2'][idx] + 1)
            except KeyError:
                quality_flags[galaxy + '_muse_ha'] = [tab['col2'][idx] + 1]

        if galaxy.upper() + '_mass' in tab['col1']:
            idx = np.where(tab['col1'] == galaxy + '_mass')[0][0]
            try:
                quality_flags[galaxy + '_muse_mass'].append(tab['col2'][idx] + 1)
            except KeyError:
                quality_flags[galaxy + '_muse_mass'] = [tab['col2'][idx] + 1]

for galaxy in alma_galaxies:
    has_alma_pattern_speed = False
    has_muse_mass_pattern_speed = False
    has_muse_ha_pattern_speed = False

    galaxy_row = galaxy_table[galaxy_table['name'] == galaxy]

    # Pull out distance and inclination from the master table, so users can convert back if necessary
    dist = galaxy_row['dist'][0]
    incl = galaxy_row['orient_incl'][0]
    pgc = galaxy_row['pgc'][0]
    pa = galaxy_row['orient_posang'][0]
    bar_r = galaxy_row['morph_bar_r'][0]
    if not np.isnan(bar_r):
        has_bar = 1
    else:
        has_bar = 0

    galaxy = galaxy.upper()
    print(galaxy)

    # Take the modal quality flag
    try:
        alma_co_q, count = mode(quality_flags[galaxy + '_alma'])

        # If everyone disagrees, take the max
        if np.any(count == 1):
            alma_co_q = np.nanmax(quality_flags[galaxy + '_alma'])
        else:
            alma_co_q = alma_co_q[0]
    except KeyError:
        alma_co_q = np.nan
    try:
        muse_mass_q, count = mode(quality_flags[galaxy + '_muse_mass'])
        if np.any(count == 1):
            alma_co_q = np.nanmax(quality_flags[galaxy + '_muse_mass'])
        else:
            muse_mass_q = muse_mass_q[0]
    except KeyError:
        muse_mass_q = np.nan
    try:
        muse_ha_q, count = mode(quality_flags[galaxy + '_muse_ha'])
        if np.any(count == 1):
            muse_ha_q = np.nanmax(quality_flags[galaxy + '_muse_ha'])
        else:
            muse_ha_q = muse_ha_q[0]
    except KeyError:
        muse_ha_q = np.nan

    # Read in the ALMA pattern speed

    try:

        file_name = alma_output + alma_version + '/'

        if not mask_outside_bar:
            file_name += 'nobmask/'
        else:
            file_name += 'bmask/'

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

        has_muse_mass_pattern_speed = True
    except OSError:
        om_muse_mass, om_muse_mass_up, om_muse_mass_down = np.nan, np.nan, np.nan

    # Read in the MUSE Halpha pattern speed

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

        has_muse_ha_pattern_speed = True
    except OSError:
        om_muse_ha, om_muse_ha_up, om_muse_ha_down = np.nan, np.nan, np.nan

    if not has_alma_pattern_speed and not has_muse_mass_pattern_speed and not has_muse_ha_pattern_speed:
        continue

    dict_keys = ['GALAXY', 'PGC', 'DIST', 'INCL', 'PA', 'HAS_BAR',
                 'OM_P_ALMA', 'OM_P_ALMA_ERR_UP', 'OM_P_ALMA_ERR_DOWN', 'OM_P_ALMA_QUAL',
                 'OM_P_MUSE_MASS', 'OM_P_MUSE_MASS_ERR_UP', 'OM_P_MUSE_MASS_ERR_DOWN', 'OM_P_MUSE_MASS_QUAL',
                 'OM_P_MUSE_HA', 'OM_P_MUSE_HA_ERR_UP', 'OM_P_MUSE_HA_ERR_DOWN', 'OM_P_MUSE_HA_QUAL',
                 ]
    dict_values = [galaxy, pgc, dist, incl, pa, has_bar,
                   om_alma, om_alma_up, om_alma_down, alma_co_q,
                   om_muse_mass, om_muse_mass_up, om_muse_mass_down, muse_mass_q,
                   om_muse_ha, om_muse_ha_up, om_muse_ha_down, muse_ha_q,
                   ]

    for i, key in enumerate(dict_keys):

        if key in data.keys():
            data[key].append(dict_values[i])
        else:
            data[key] = [dict_values[i]]

    if not keys:
        keys = list(data.keys())

    # Check for any ALMA radii

    data, alma_corot_keys = add_radii(resonance_output + galaxy + '_alma_corot.txt', data, alma_corot_keys)
    data, alma_ilr_keys = add_radii(resonance_output + galaxy + '_alma_ilr.txt', data, alma_ilr_keys,
                                    radii_name='R_ILR_ALMA')
    data, alma_olr_keys = add_radii(resonance_output + galaxy + '_alma_olr.txt', data, alma_olr_keys,
                                    radii_name='R_OLR_ALMA')

    # Check for any MUSE mass radii

    data, muse_mass_corot_keys = add_radii(resonance_output + galaxy + '_muse_mass_corot.txt', data,
                                           muse_mass_corot_keys, radii_name='R_CR_MUSE_MASS')
    data, muse_mass_ilr_keys = add_radii(resonance_output + galaxy + '_muse_mass_ilr.txt', data, muse_mass_ilr_keys,
                                         radii_name='R_ILR_MUSE_MASS')
    data, muse_mass_olr_keys = add_radii(resonance_output + galaxy + '_muse_mass_olr.txt', data, muse_mass_olr_keys,
                                         radii_name='R_OLR_MUSE_MASS')

    # Check for any MUSE Ha radii

    data, muse_ha_corot_keys = add_radii(resonance_output + galaxy + '_muse_ha_corot.txt', data, muse_ha_corot_keys,
                                         radii_name='R_CR_MUSE_HA')
    data, muse_ha_ilr_keys = add_radii(resonance_output + galaxy + '_muse_ha_ilr.txt', data, muse_ha_ilr_keys,
                                       radii_name='R_ILR_MUSE_HA')
    data, muse_ha_olr_keys = add_radii(resonance_output + galaxy + '_muse_ha_olr.txt', data, muse_ha_olr_keys,
                                       radii_name='R_OLR_MUSE_HA')

for key_to_extend in [alma_corot_keys, alma_ilr_keys, alma_olr_keys,
                      muse_mass_corot_keys, muse_mass_ilr_keys, muse_mass_olr_keys,
                      muse_ha_corot_keys, muse_ha_ilr_keys, muse_ha_olr_keys]:
    keys.extend(key_to_extend)

table = Table()

for key in keys:

    column = data[key]

    if 'DIST' in key:
        column = column * u.Mpc
    elif 'INCL' in key:
        column = column * u.degree
    elif 'PA' in key:
        column = column * u.degree
    elif 'OM_P' in key and 'QUAL' not in key:
        column = column * u.km / u.s / u.kpc
    elif 'R_CR' in key or 'R_ILR' in key or 'R_OLR' in key:
        column = column * u.kpc

    table.add_column(column, name=key)

table.write(output_folder + 'pattern_speed_table_' + pattern_speed_version + '.fits',
            format='fits', overwrite=True)

# Write out a table for just the good stellar masses

names = ['Galaxy', 'PGC', 'D', '$i$', 'PA', 'Bar?']
units = {'D': u.Mpc,
         '$i$': u.deg,
         'PA': u.deg,
         r'$\Omega_{\rm P}$': u.km * u.s ** -1 * u.kpc ** -1,
         r'$R_{\rm CR}$': u.kpc
         }
col_align = 'c' * len(names)

rows = table[(table['OM_P_MUSE_MASS_QUAL'] == 1) | (table['OM_P_MUSE_MASS_QUAL'] == 2)]

latex_table = Table()
latex_table.add_column(rows['GALAXY'], name='Galaxy')
latex_table.add_column(rows['PGC'], name='PGC')
latex_table.add_column(rows['DIST'], name='D')
latex_table.add_column(rows['INCL'], name='i')
latex_table.add_column(rows['PA'], name='PA')
latex_table.add_column(rows['HAS_BAR'], name='Bar?')

latex_muse_mass_om_p, names = texify_om_p(latex_table, names, rows, data_name='MUSE_MASS', short_name='MM')
col_align += 'cc'

latex_table, col_align, names = texify_resonances(muse_mass_corot_keys[:1], rows, latex_table, col_align,
                                                  names, r'$R_{\rm CR%d, MM}$', n_cols=n_corot_cols)
names[-3] = '$\\Omega_{\\rm P}$'
names[-2] = 'Q'
names[-1] = '$R_{\\rm CR}$'

for key in units.keys():
    units[key] = units[key].to_string('latex')

n_galaxies = len(latex_table)

ascii.write(latex_table,
            output_folder + 'pattern_speed_table_muse_mass_only.tex',
            col_align=col_align,
            # sys.stdout,
            names=names,
            fill_values=[('nan', r'\nodata'),
                         ('-9223372036854775808', r'\nodata')],
            format='aastex',
            overwrite=True,
            latexdict={'units': units,
                       'tabletype': 'deluxetable*',
                       'preamble': r'\tablecaption{Pattern speeds and co-rotation radii for the ten well-constrained'
                                   r' stellar mass pattern speeds.}'
                                   '\n\label{table:well_constrained_pattern_speeds}'},
            )

# TeXiFy this table for public consumption

names = ['Galaxy', 'PGC', 'D', '$i$', 'PA', 'Bar?']
units = {'D': u.Mpc,
         '$i$': u.deg,
         'PA': u.deg,
         r'$\Omega_{\rm clump, A}$': u.km * u.s ** -1 * u.kpc ** -1,
         r'$\Omega_{\rm P, MM}$': u.km * u.s ** -1 * u.kpc ** -1,
         r'$\Omega_{\rm clump, MH{\alpha}}$': u.km * u.s ** -1 * u.kpc ** -1,
         }

col_align = 'c' * len(names)

latex_table = Table()
latex_table.add_column(table['GALAXY'], name='Galaxy')
latex_table.add_column(table['PGC'], name='PGC')
latex_table.add_column(table['DIST'], name='D')
latex_table.add_column(table['INCL'], name='i')
latex_table.add_column(table['PA'], name='PA')
latex_table.add_column(table['HAS_BAR'], name='Bar?')

# Include ALMA stuff

latex_table, names = texify_om_p(latex_table, names, table, data_name='ALMA', short_name='A')
col_align += 'cc'

latex_table, col_align, names = texify_resonances(alma_corot_keys, table, latex_table, col_align,
                                                  names, r'$R_{\rm CR%d, A}$', n_cols=n_corot_cols)
latex_table, col_align, names = texify_resonances(alma_ilr_keys, table, latex_table, col_align,
                                                  names, r'$R_{\rm ILR%d, A}$', n_cols=0)
latex_table, col_align, names = texify_resonances(alma_olr_keys, table, latex_table, col_align,
                                                  names, r'$R_{\rm OLR%d, A}$', n_cols=0)

# Include the MUSE mass stuff

latex_muse_mass_om_p, names = texify_om_p(latex_table, names, table, data_name='MUSE_MASS', short_name='MM')
col_align += 'cc'

latex_table, col_align, names = texify_resonances(muse_mass_corot_keys, table, latex_table, col_align,
                                                  names, r'$R_{\rm CR%d, MM}$', n_cols=n_corot_cols)
latex_table, col_align, names = texify_resonances(muse_mass_ilr_keys, table, latex_table, col_align,
                                                  names, r'$R_{\rm ILR%d, MM}$', n_cols=0)
latex_table, col_align, names = texify_resonances(muse_mass_olr_keys, table, latex_table, col_align,
                                                  names, r'$R_{\rm OLR%d, MM}$', n_cols=0)

# Include MUSE Halpha stuff

latex_table, names = texify_om_p(latex_table, names, table, data_name='MUSE_HA', short_name=r'MH{\alpha}')
col_align += 'cc'

latex_table, col_align, names = texify_resonances(muse_ha_corot_keys, table, latex_table, col_align,
                                                  names, r'$R_{\rm CR%d, MH{\alpha}}$', n_cols=n_corot_cols)
latex_table, col_align, names = texify_resonances(muse_ha_ilr_keys, table, latex_table, col_align,
                                                  names, r'$R_{\rm ILR%d, MH{\alpha}}$', n_cols=0)
latex_table, col_align, names = texify_resonances(muse_ha_olr_keys, table, latex_table, col_align,
                                                  names, r'$R_{\rm OLR%d, MH{\alpha}}$', n_cols=0)


for key in units.keys():
    units[key] = units[key].to_string('latex')

n_galaxies = len(latex_table)

ascii.write(latex_table,
            output_folder + 'pattern_speed_table_' + pattern_speed_version + '.tex',
            col_align=col_align,
            # sys.stdout,
            names=names,
            fill_values=[('nan', r'\nodata'),
                         ('-9223372036854775808', r'\nodata')],
            format='aastex',
            overwrite=True,
            latexdict={'units': units,
                       'tabletype': 'deluxetable*',
                       'preamble': r"\tablecaption{``Pattern speeds'' ($\Omega_p$ for stellar mass, "
                                   r'$\Omega_{\rm clump}$ for ISM tracers) and inferred resonance locations for %d '
                                   r'galaxies. A description of the column names is given in the Table notes.}'
                                   '\n\label{table:pattern_speeds}' % n_galaxies},
            )

# In each case, we use a subscript to refer to the tracer in\\ question --'
#                                    r'MM for MUSE stellar mass, MH$\alpha$ for MUSE-H$\alpha$ measurements, and A for '
#                                    r'ALMA CO measurements. We also show the quality flags (indicated as Q), and \\'
#                                    r'co-rotation radii ($R_{\rm CR}$), and outer/inner Lindblad radii '
#                                    r'($R_{\rm O/ILR}$). For brevity, only the first co-rotation radius is shown, and '
#                                    r'Lindblad radii are hidden.

# Finally, we need to swap 'colhead' to 'nocolhead' for the hidden columns.

f = open(output_folder + 'pattern_speed_table_' + pattern_speed_version + '.tex', 'r')

col_align_edit = col_align + col_align

lines = f.readlines()
line = lines[3]

new_line = ''
columns = line.split('colhead')

for i, string in enumerate(columns[:-1]):

    new_line += string
    if col_align_edit[i] == 'c':
        new_line += 'colhead'
    else:
        new_line += 'nocolhead'

new_line += columns[-1]

lines[3] = new_line

# Add in the table notes

lines.insert(-1,
             r'\tablecomments{For each galaxy, we list its NGC number, PGC number, distance, inclination, position '
             r'angle, and whether it has a bar (1 for yes, 0 for no). We use a subscript to refer to the tracer in '
             r'question -- MM for MUSE stellar mass, '
             r'MH$\alpha$ for MUSE-H$\alpha$ measurements, and A for ALMA CO measurements. We also show the quality '
             r'flags (indicated as Q), and co-rotation radii ($R_{\rm CR}$), and outer/inner Lindblad radii '
             r'($R_{\rm O/ILR}$). For brevity, only the first co-rotation radius is shown, and Lindblad radii are '
             r'hidden.}')

f = open(output_folder + 'pattern_speed_table_' + pattern_speed_version + '.tex', 'w')
f.writelines(lines)

f.close()

print('Complete!')
