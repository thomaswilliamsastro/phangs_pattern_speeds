# -*- coding: utf-8 -*-
"""
Calculate correlations between various parameters from this study and galaxy parameters

@author: Tom Williams
"""

import os
import time

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from astropy.table import Table
from matplotlib import gridspec
from scipy.stats import kendalltau

from vars import phangs_folder, output_folder, alma_galaxies, pattern_speed_version, galaxy_table, plot_folder


def scatterplot_with_correlation(x_data, y_data, x_lims=None, y_lims=None, c_lims=None,
                                 different_symbols=None, edge_colours=None,
                                 x_label=None, show_x_label=True, y_label=None, show_y_label=True,
                                 subplot=(1, 1, 1),
                                 c='k', c_label=None, use_errors=False, x_data_errs=None, y_data_errs=None):
    """Create scatterplot with Kendall Tau/p-value text.
    """

    x_data = np.array(x_data)
    y_data = np.array(y_data)

    ax = plt.subplot(subplot)

    if not edge_colours:
        edge_colours = ['none'] * len(x_data)
    edge_colours = np.array(edge_colours)

    # Filter any NaNs

    nan_idx = np.where((np.isnan(x_data) == False) & (np.isnan(y_data) == False))

    n = len(nan_idx[0])

    if not use_errors:
        tau = kendalltau(np.asarray(x_data)[nan_idx], np.asarray(y_data)[nan_idx])
        tau = [tau[0], tau[1]]
        if tau[1] < 0.01:
            tau[1] = '<0.01'
        else:
            tau[1] = '=%.2f' % tau[1]
    else:
        n_draws = 10000
        tau_mc = np.zeros(n_draws)
        for i in range(n_draws):
            x_data_mc = x_data + np.random.normal(scale=x_data_errs)
            y_data_mc = y_data + np.random.normal(scale=y_data_errs)
            tau = kendalltau(np.asarray(x_data_mc)[nan_idx], np.asarray(y_data_mc)[nan_idx])
            tau_mc[i] = tau[0]

        tau_mc_median = np.nanmedian(tau_mc)
        tau_mc_err = np.nanpercentile(tau_mc, 84) - tau_mc_median

    if type(c) is not str:
        c_nan_idx = np.where(np.isnan(c) == True)

        if len(c_nan_idx[0] > 0):
            plt.scatter(np.asarray(x_data)[c_nan_idx], np.asarray(y_data)[c_nan_idx],
                        edgecolors='k', facecolor='none', alpha=0.5)

    if c_lims is None:
        if different_symbols is None:
            scatter = ax.scatter(x_data, y_data, c=c, edgecolors=edge_colours)
        else:
            mask_nan = np.where(np.isnan(different_symbols))
            mask_non_nan = np.where(~np.isnan(different_symbols))
            scatter = ax.scatter(x_data[mask_nan], y_data[mask_nan], marker='s', c=c, edgecolors=edge_colours[mask_nan])
            ax.scatter(x_data[mask_non_nan], y_data[mask_non_nan], marker='o', c=c,
                       edgecolors=edge_colours[mask_non_nan])
    else:
        c = np.array(c)
        if different_symbols is None:
            scatter = ax.scatter(x_data, y_data, c=c, vmin=c_lims[0], vmax=c_lims[1], edgecolors=edge_colours)
        else:
            mask_nan = np.where(np.isnan(different_symbols))
            mask_non_nan = np.where(~np.isnan(different_symbols))
            scatter = ax.scatter(x_data[mask_nan], y_data[mask_nan], marker='s', c=c[mask_nan], vmin=c_lims[0],
                                 vmax=c_lims[1], edgecolors=edge_colours[mask_nan])
            ax.scatter(x_data[mask_non_nan], y_data[mask_non_nan],
                       marker='o', c=c[mask_non_nan], vmin=c_lims[0], vmax=c_lims[1],
                       edgecolors=edge_colours[mask_non_nan])

    if not use_errors:
        ax.text(0.95, 0.95,
                '$N= %d$,\n$\\tau=%.2f, p%s$' % (n, tau[0], tau[1]),
                ha='right', va='top', transform=ax.transAxes,
                bbox=dict(facecolor='white', alpha=0.5))
    else:
        ax.text(0.95, 0.95,
                '$N= %d$,\n$\\tau=%.2f\pm%.2f$' % (n, tau_mc_median, tau_mc_err),
                ha='right', va='top', transform=ax.transAxes,
                bbox=dict(facecolor='white', alpha=0.5))

    if x_label is not None:
        plt.xlabel(x_label)
    if y_label is not None:
        plt.ylabel(y_label)

    if not show_x_label:
        ax.tick_params(labelleft=False)
    if not show_y_label:
        ax.tick_params(labelbottom=False)

    if x_lims is not None:
        plt.xlim(x_lims)
    if y_lims is not None:
        plt.ylim(y_lims)

    if c_label is not None:
        plt.colorbar(scatter, label=c_label)

    return ax, scatter


def multi_panel_plot(x_quantities, y_quantities, parameters, limits, labels,
                     c_quantities='alma_incl', show_c_left=False, different_symbols=None, edge_colours=None,
                     cax_loc=None, background_cloud=None, use_errors=False):
    if cax_loc is None:
        cax_loc = [0.9125, 0.1, 0.0125, 0.8]
    fig = plt.figure(figsize=(3 * len(x_quantities), 2 * len(y_quantities)))

    fig.subplots_adjust(left=0.1, right=0.8, bottom=0.1, top=0.9)

    gs1 = gridspec.GridSpec(len(y_quantities), len(x_quantities))
    gs1.update(hspace=0, wspace=0)

    subplot_n = 0

    for y_quantity in y_quantities:

        for x_quantity in x_quantities:

            if not c_quantities:

                c = 'k'
                c_lims = None
                c_label = None

            else:

                if type(c_quantities) == dict:

                    c_quantity = c_quantities[x_quantity]

                else:
                    c_quantity = c_quantities

                c = parameters[c_quantity]
                c_lims = limits[c_quantity]
                c_label = labels[c_quantity]

            if not y_quantity == y_quantities[-1]:
                ax, _ = scatterplot_with_correlation(parameters[x_quantity], parameters[y_quantity],
                                                     x_lims=limits[x_quantity], y_lims=limits[y_quantity],
                                                     c=c, c_lims=c_lims,
                                                     y_label=labels[y_quantity],
                                                     different_symbols=different_symbols, edge_colours=edge_colours,
                                                     subplot=gs1[subplot_n],
                                                     use_errors=use_errors,
                                                     x_data_errs=parameters[x_quantity + '_err'],
                                                     y_data_errs=parameters[y_quantity + '_err'])

            else:
                ax, scatter = scatterplot_with_correlation(parameters[x_quantity], parameters[y_quantity],
                                                           x_lims=limits[x_quantity], y_lims=limits[y_quantity],
                                                           c=c, c_lims=c_lims,
                                                           x_label=labels[x_quantity],
                                                           y_label=labels[y_quantity],
                                                           different_symbols=different_symbols,
                                                           edge_colours=edge_colours,
                                                           subplot=gs1[subplot_n],
                                                           use_errors=use_errors,
                                                           x_data_errs=parameters[x_quantity + '_err'],
                                                           y_data_errs=parameters[y_quantity + '_err'])

                if x_quantity == x_quantities[0] and show_c_left:
                    cax = plt.axes([-0.05, 0.1, 0.0125, 0.8])

                    cbar_ax = plt.colorbar(scatter, cax=cax, label=c_label)

                    cax.yaxis.set_ticks_position('left')
                    cax.yaxis.set_label_position('left')

                if x_quantity == x_quantities[-1]:
                    cax = plt.axes(cax_loc)

                    cbar_ax = plt.colorbar(scatter, cax=cax, label=c_label)

            if x_quantity == x_quantities[-1]:
                ax.yaxis.set_label_position("right")
                ax.yaxis.tick_right()

            if x_quantity not in [x_quantities[0], x_quantities[-1]]:
                ax.yaxis.set_ticklabels([])
                plt.ylabel('')

            if background_cloud is not None:
                ax.scatter(parameters[x_quantity.replace(background_cloud[0], background_cloud[1])],
                           parameters[y_quantity.replace(background_cloud[0], background_cloud[1])],
                           c='k', marker='x', alpha=0.25, zorder=-1)

            subplot_n += 1

    return fig


def corner(quantities, parameters, limits, labels, different_symbols=None, edge_colours=None,
           use_errors=False):
    fig = plt.figure(figsize=(3 * (len(quantities) - 1), 2 * (len(quantities) - 1)))

    fig.subplots_adjust(left=0.2, right=0.8, bottom=0.1, top=0.9)

    gs1 = gridspec.GridSpec(len(quantities) - 1, len(quantities) - 1)
    gs1.update(hspace=0, wspace=0)

    for i, x_quantity in enumerate(quantities):

        subplot_n = (len(quantities) - 1) * i + i

        for j, y_quantity in enumerate(quantities[i + 1:]):

            already_plotted = False

            if subplot_n > (len(quantities) - 1) ** 2 - len(quantities):

                ax, _ = scatterplot_with_correlation(parameters[x_quantity], parameters[y_quantity],
                                                     x_lims=limits[x_quantity], y_lims=limits[y_quantity],
                                                     x_label=labels[x_quantity],
                                                     subplot=gs1[subplot_n],
                                                     different_symbols=different_symbols, edge_colours=edge_colours,
                                                     use_errors=use_errors,
                                                     x_data_errs=parameters[x_quantity + '_err'],
                                                     y_data_errs=parameters[y_quantity + '_err']
                                                     )

                if subplot_n % (len(quantities) - 1) != 0:
                    ax.yaxis.set_ticklabels([])
                    plt.ylabel('')

                already_plotted = True

            if subplot_n % (len(quantities) - 1) == 0:
                ax, _ = scatterplot_with_correlation(parameters[x_quantity], parameters[y_quantity],
                                                     x_lims=limits[x_quantity], y_lims=limits[y_quantity],
                                                     y_label=labels[y_quantity],
                                                     subplot=gs1[subplot_n],
                                                     different_symbols=different_symbols, edge_colours=edge_colours,
                                                     use_errors=use_errors,
                                                     x_data_errs=parameters[x_quantity + '_err'],
                                                     y_data_errs=parameters[y_quantity + '_err']
                                                     )

                already_plotted = True

            if not already_plotted:
                ax, _ = scatterplot_with_correlation(parameters[x_quantity], parameters[y_quantity],
                                                     x_lims=limits[x_quantity], y_lims=limits[y_quantity],
                                                     subplot=gs1[subplot_n],
                                                     different_symbols=different_symbols, edge_colours=edge_colours,
                                                     use_errors=use_errors,
                                                     x_data_errs=parameters[x_quantity + '_err'],
                                                     y_data_errs=parameters[y_quantity + '_err']
                                                     )

                ax.yaxis.set_ticklabels([])
                ax.xaxis.set_ticklabels([])
                plt.ylabel('')

            subplot_n += len(quantities) - 1

    return fig


matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 14

os.chdir(phangs_folder)

tracers_to_use = ['alma', 'muse_mass']

start = time.time()

pattern_speeds_t = Table.read(output_folder + 'pattern_speed_table_' + pattern_speed_version + '.fits')

# Pull in scale lengths

s4g_galaxy, scale_lengths, b_t_ratios = np.loadtxt('lit_tables/s4g_scale_lengths.txt',
                                                   unpack=True, dtype=str)

# And bar strengths

b_s_galaxy, bar_strengths = np.loadtxt('lit_tables/bar_strengths_diazgarcia16.dat', unpack=True, dtype=str,
                                       usecols=(0, 2))

# The table containing aperture corrected CO luminosities

co_lum_table = Table.read('documents/galaxy_table.fits')

err_cutoff = 0.5

parameters = {}
quantities = ['delta_pa', 'incl', 'morph', 'mstar', 'hr', 'r25', 'b_t',
              'r_bar', 'r_cr', 'r', 'r_err', 'om_p',
              's_bar', 'f_mol', 'vr_inf', 'r_t', 'spiral_arms',
              'galaxies', 'q_flags']

old_sample_table = Table.read('documents/phangs_sample_table_v1p4.fits')

for quantity in quantities:
    for tracer_to_use in tracers_to_use:
        parameters[tracer_to_use + '_' + quantity] = []
        parameters[tracer_to_use + '_' + quantity + '_err'] = []

for tracer_to_use in tracers_to_use:
    for galaxy in alma_galaxies:

        # From the master table, pull out various parameters

        gal_row = galaxy_table[galaxy_table['name'] == galaxy]
        co_lum_row = co_lum_table[co_lum_table['NAME'] == galaxy.ljust(10)]

        dist = gal_row

        gal_pa = gal_row['orient_posang'][0]
        gal_pa_err = gal_row['orient_posang_unc'][0]

        morph = gal_row['morph_t'][0]
        morph_err = gal_row['morph_t_unc'][0]
        m_star = gal_row['props_mstar'][0]
        m_star_err = 10 ** (np.log10(m_star) + gal_row['props_mstar_unc'][0]) - m_star
        r_25 = gal_row['size_r25'][0]
        r_25_err = 0.1 * r_25
        v_rot = gal_row['rotcur_v0'][0]
        v_rot_err = gal_row['rotcur_v0_unc'][0]
        r_t = gal_row['rotcur_rt'][0]
        r_t_err = gal_row['rotcur_rt_unc'][0]

        spiral_arms = gal_row['morph_spiral_arms'][0] + gal_row['morph_multi_arm'][0]
        if spiral_arms >= 1:
            spiral_arms = 1
        else:
            spiral_arms = np.nan

        bar_r = gal_row['morph_bar_r'][0]
        bar_pa = gal_row['morph_bar_pa'][0]

        # Calculate the molecular gas mass and fraction from depletion timescale, SFR, and Mstar

        try:
            m_mol = co_lum_row['TDEP_PHANGS'][0] * co_lum_row['SFR_Z0MGS'][0]
        except IndexError:
            m_mol = np.nan

        f_mol = m_mol / (m_mol + m_star)
        f_mol_err = 0.1 * f_mol

        # Pull out galaxy row from pattern speeds

        row = pattern_speeds_t[pattern_speeds_t['GALAXY'] == galaxy.upper()]
        if len(row) == 0:
            continue

        # Convert bar to physical distance

        dist = row['DIST'][0]
        dist_err = 10 ** (np.log10(dist) + gal_row['dist_unc'][0]) - dist

        bar_r *= np.pi / 180 * 1 / 3600 * dist * 1e3

        incl = row['INCL'][0]
        incl_err = gal_row['orient_incl_unc'][0]

        # Account for distance differences

        old_sample_table_idx = np.where(old_sample_table['NAME'] == galaxy.upper().ljust(10, ' '))[0][0]
        old_dist = old_sample_table['DIST'][old_sample_table_idx]

        r_t *= dist / old_dist
        r_t_err *= dist / old_dist

        # Deproject bar length

        delta_pa = bar_pa - gal_pa

        bar_r = bar_r * np.sqrt((np.cos(np.radians(delta_pa)) ** 2 +
                                 np.sin(np.radians(delta_pa)) ** 2 / np.cos(np.radians(incl)) ** 2))
        bar_r_err = 0.1 * bar_r

        delta_pa = np.abs(delta_pa)

        if delta_pa > 180:
            delta_pa -= 180
        if delta_pa > 90:
            delta_pa = 180 - delta_pa

        galaxy = galaxy.upper()

        # Pull in scale length and convert to physical distance

        s4g_idx = np.where(s4g_galaxy == galaxy)
        scale_length = float(scale_lengths[s4g_idx[0][0]])

        scale_length *= np.pi / 180 * 1 / 3600 * dist * 1e3
        scale_length_err = scale_length * dist_err / dist

        # Bulge-to-disk ratio

        b_t = float(b_t_ratios[s4g_idx[0][0]])
        b_t_err = 0.1 * b_t

        # Convert r_25 to physical distance

        r_25 *= np.pi / 180 * 1 / 3600 * dist * 1e3

        # Bar strength

        b_s_idx = np.where(b_s_galaxy == galaxy)[0]
        if len(b_s_idx) == 0:
            b_s = np.nan
        else:
            b_s = float(bar_strengths[b_s_idx])

        b_s_err = 0.1 * b_s

        # Pull the corotation radii and quality flags from the pattern speeds table

        om_p = row['OM_P_' + tracer_to_use.upper()][0]
        om_p_err = np.mean([row['OM_P_' + tracer_to_use.upper() + '_ERR_UP'][0],
                            row['OM_P_' + tracer_to_use.upper() + '_ERR_DOWN'][0]])
        q = row['OM_P_' + tracer_to_use.upper() + '_QUAL'][0]

        # Only use quality 1/2 pattern speeds

        if q in [1, 2]:

            # if tracer_to_use == 'muse_mass':
            #     print(r_25)

            om_p_err_percent = om_p_err / om_p

            corot = np.nan
            corot_err = np.nan

            if tracer_to_use == 'alma':
                corot_range = 4
            else:
                corot_range = 2

            for i in range(corot_range):

                if row['R_CR_' + tracer_to_use.upper() + '_' + str(i + 1)][0] > 1:
                    corot = row['R_CR_' + tracer_to_use.upper() + '_' + str(i + 1)][0]
                    corot_err = row['R_CR_' + tracer_to_use.upper() + '_' + str(i + 1) + '_ERR'][0]
                    break

            if om_p_err_percent < err_cutoff:
                parameters[tracer_to_use + '_delta_pa'].append(delta_pa)
                parameters[tracer_to_use + '_delta_pa_err'].append(gal_pa_err)
                parameters[tracer_to_use + '_incl'].append(incl)
                parameters[tracer_to_use + '_incl_err'].append(incl_err)

                parameters[tracer_to_use + '_morph'].append(morph)
                parameters[tracer_to_use + '_morph_err'].append(morph_err)
                parameters[tracer_to_use + '_mstar'].append(np.log10(m_star) - np.log10(r_25 ** 2))
                parameters[tracer_to_use + '_mstar_err'].append(np.log10(m_star_err) - np.log10(r_25 ** 2))
                parameters[tracer_to_use + '_hr'].append(scale_length)
                parameters[tracer_to_use + '_hr_err'].append(scale_length_err)
                parameters[tracer_to_use + '_r25'].append(r_25)
                parameters[tracer_to_use + '_r25_err'].append(r_25_err)
                parameters[tracer_to_use + '_b_t'].append(b_t)
                parameters[tracer_to_use + '_b_t_err'].append(b_t_err)

                parameters[tracer_to_use + '_r_bar'].append(bar_r / r_25)
                parameters[tracer_to_use + '_r_bar_err'].append(bar_r_err)
                parameters[tracer_to_use + '_r_cr'].append(corot)
                parameters[tracer_to_use + '_r_cr_err'].append(corot_err)
                parameters[tracer_to_use + '_r'].append(corot / bar_r)
                parameters[tracer_to_use + '_r_err'].append((corot / bar_r) *
                                                            np.sqrt(
                                                                (corot_err / corot) ** 2 + (0.2 * bar_r / bar_r) ** 2))
                parameters[tracer_to_use + '_om_p'].append(om_p)
                parameters[tracer_to_use + '_om_p_err'].append(om_p_err)
                parameters[tracer_to_use + '_s_bar'].append(b_s)
                parameters[tracer_to_use + '_s_bar_err'].append(b_s_err)

                parameters[tracer_to_use + '_f_mol'].append(f_mol)
                parameters[tracer_to_use + '_f_mol_err'].append(f_mol_err)
                parameters[tracer_to_use + '_vr_inf'].append(v_rot)
                parameters[tracer_to_use + '_vr_inf_err'].append(v_rot_err)
                parameters[tracer_to_use + '_r_t'].append(r_t)
                parameters[tracer_to_use + '_r_t_err'].append(r_t_err)
                parameters[tracer_to_use + '_spiral_arms'].append(spiral_arms)
                parameters[tracer_to_use + '_q_flags'].append(q)
                parameters[tracer_to_use + '_galaxies'].append(galaxy)

            else:
                print('%s fails err cutoff' % galaxy)

main_tracer = 'muse_mass'
background_tracer = 'alma'

# pattern speed found for bars, spirals (without bars), spirals (with bars), no bar and no spiral
# pattern speed not found for bars, spirals (without bars), spirals (with bars), no bar and no spiral

# Number statistics
# mask = np.where( (np.asarray(parameters['alma_q_flags']) == 3) & (np.asarray(parameters['alma_r_bar']) > 0) )
# print(len(mask[0]))
# no
#
# mask = np.where( (np.asarray(parameters['alma_spiral_arms']) == 1) & (np.asarray(parameters['alma_r_bar']) > 0) )
# print(np.asarray(parameters['alma_q_flags'])[mask])

om_p = np.array(parameters[main_tracer + '_om_p'])
m_star = np.array(parameters[main_tracer + '_mstar'])
r25 = np.array(parameters[main_tracer + '_r25'])
m_star += np.log10(r25 ** 2)
m_star = 10 ** m_star

# for i, alma_r in enumerate(parameters['alma_r']):
#
#     alma_r_err = parameters['alma_r_err'][i]
#
#     if not alma_r - alma_r_err <= 1.2 <= alma_r + alma_r_err and not np.isnan(alma_r):
#
#         print(alma_r)

median_om_p_m_star = np.nanmedian(om_p / m_star)
om_p_m_star_err_up = np.nanpercentile(om_p / m_star, 84) - median_om_p_m_star
om_p_m_star_err_down = median_om_p_m_star - np.nanpercentile(om_p / m_star, 16)
print('Mass-weighted average pattern speed: $%.2e^{+%.2e}_{-%.2e}$' %
      (median_om_p_m_star, om_p_m_star_err_up, om_p_m_star_err_down))

c_incl = np.nanpercentile(np.asarray(parameters[main_tracer + '_incl']), [1, 99])
c_delta_pa = np.nanpercentile(np.asarray(parameters[main_tracer + '_delta_pa']), [1, 99])

parameter_to_edgecolor = main_tracer + '_spiral_arms'

edge_colours = []
for parameter in parameters[parameter_to_edgecolor]:

    if parameter > 0:
        edge_colours.append('r')
    else:
        edge_colours.append('none')

# edge_colours = []
# for i, parameter in enumerate(parameters['alma_spiral_arms']):
#
#     if np.isnan(parameter) and np.isnan(parameters['alma_r_bar'][i]):
#         edge_colours.append('r')
#     else:
#         edge_colours.append('none')

# PLOT ONE: Distribution of fancy R for whole galaxy sample

print(np.nanpercentile(np.array(parameters['muse_mass_om_p_err']) / np.array(parameters['muse_mass_om_p']),
                       [16, 50, 84]))

print(np.nanpercentile(parameters[main_tracer + '_r'], [16, 50, 84]))

plt.figure(figsize=(4, 3))

sns.kdeplot(parameters[main_tracer + '_r'], bw='silverman', color='k', shade=True,
            clip=[np.nanmin(parameters[main_tracer + '_r']), np.nanmax(parameters[main_tracer + '_r'])])

plt.axvline(np.nanmedian(parameters[main_tracer + '_r']), c='k', ls='-')
plt.axvline(np.nanpercentile(parameters[main_tracer + '_r'], 16), c='k', ls='--')
plt.axvline(np.nanpercentile(parameters[main_tracer + '_r'], 84), c='k', ls='--')
plt.axvline(1.2, c='r', ls='-')

xlims = plt.xlim()

plt.xlim([0, xlims[-1]])

plt.xlabel(r'$\mathcal{R}$')
plt.ylabel('Probability Density')

# plt.show()

plt.savefig(plot_folder + 'bar_corot_r_ratio_' + main_tracer + '.png',
            bbox_inches='tight')
plt.savefig(plot_folder + 'bar_corot_r_ratio_' + main_tracer + '.pdf',
            bbox_inches='tight')

plt.close()

ratio = 2 * np.array(parameters[main_tracer + '_r_cr']) / np.array(parameters[main_tracer + '_hr'])

print(np.nanpercentile(ratio, [16, 50, 84]))

plt.figure(figsize=(4, 3))

sns.kdeplot(ratio, bw='silverman', color='k', shade=True)
plt.axvline(np.nanmedian(ratio), c='k', ls='-')
plt.axvline(np.nanpercentile(ratio, 16), c='k', ls='--')
plt.axvline(np.nanpercentile(ratio, 84), c='k', ls='--')
plt.axvline(3, c='r', ls='-')

xlims = plt.xlim()

plt.xlim([0, xlims[-1]])

plt.xlabel(r'$R_\mathrm{CR}/h_R$')
plt.ylabel('Probability Density')

# plt.show()

plt.savefig(plot_folder + 'r_cr_scale_length_' + main_tracer + '.png',
            bbox_inches='tight')
plt.savefig(plot_folder + 'r_cr_scale_length_' + main_tracer + '.pdf',
            bbox_inches='tight')

plt.close()

limits = {main_tracer + '_delta_pa': [-5, 95],
          main_tracer + '_incl': [5, 85],
          main_tracer + '_morph': [0.25, 7.75],
          main_tracer + '_hr': [1.25, 8.75],
          main_tracer + '_om_p': [15, 85],
          main_tracer + '_r_cr': [-0.5, 11.5],
          # 'alma_r_bar': [-0.5, 12.25],
          main_tracer + '_r_bar': [-0.05, 0.7],
          main_tracer + '_r': [-0.25, 5.5],
          main_tracer + '_mstar': [7.4, 8.9],
          main_tracer + '_b_t': [-0.1, 0.5],
          main_tracer + '_s_bar': [0.2, 1.0],
          main_tracer + '_f_mol': [-0.05, 0.3],
          main_tracer + '_vr_inf': [110, 360],
          main_tracer + '_r_t': [-0.1, 3.1]}

labels = {main_tracer + '_delta_pa': r'$|\Delta \mathregular{PA}|$ (deg)',
          main_tracer + '_incl': r'$i$ (deg)',
          main_tracer + '_morph': r'$T_\mathrm{Hubble}$',
          main_tracer + '_hr': r'$h_r$ (kpc)',
          main_tracer + '_om_p': r'$\Omega_\mathrm{P}~(\mathrm{km\,s^{-1}\,kpc^{-1}})$',
          main_tracer + '_r_cr': r'$R_\mathrm{CR}$ (kpc)',
          # 'alma_r_bar': r'$R_\mathrm{bar}$ (kpc)',
          main_tracer + '_r_bar': r'$R_\mathrm{bar}/R_{25}$',
          main_tracer + '_r': r'$\mathcal{R}$',
          main_tracer + '_mstar': r'$\log_{10}(M_\ast\,R_{25}^{-2} [M_\odot\,\mathrm{kpc}^{-2}])$',
          main_tracer + '_b_t': 'B/T',
          main_tracer + '_s_bar': r'$S_\mathrm{bar}$',
          main_tracer + '_f_mol': r'$f_\mathrm{mol}$',
          main_tracer + '_vr_inf': r'$v_\mathrm{r, inf}$ (km s$^{-1}$)',
          main_tracer + '_r_t': r'$R_t$ (kpc)'}

# PLOT TWO: Scatter plot of pattern speed, corotation radius, bar length, and fancy R as function of both delta PA and
# inclination. Clearly, can only define delta PA with the barred galaxies.

x_quantities = [main_tracer + '_delta_pa', main_tracer + '_incl']
y_quantities = [main_tracer + '_om_p', main_tracer + '_r_cr', main_tracer + '_r_bar', main_tracer + '_r',
                main_tracer + '_s_bar']
c_quantities = {main_tracer + '_delta_pa': main_tracer + '_incl',
                main_tracer + '_incl': main_tracer + '_delta_pa'}

fig = multi_panel_plot(x_quantities, y_quantities,
                       parameters, limits, labels, c_quantities=c_quantities, show_c_left=True, use_errors=True,
                       background_cloud=['muse_mass', 'alma'])

# plt.show()

plt.savefig(plot_folder + 'correlations_projection_' + main_tracer + '.png',
            bbox_inches='tight')
plt.savefig(plot_folder + 'correlations_projection_' + main_tracer + '.pdf',
            bbox_inches='tight')

plt.close()

# PLOT THREE: Scatter plot of pattern speed, corotation radius, bar length, and fancy R as function of galaxy
# parameters -- morphology, bulge to disk ratio, scale length, stellar mass, asymptotic velocity and molecular gas
# fraction

x_quantities = [main_tracer + '_morph', main_tracer + '_b_t', main_tracer + '_hr', main_tracer + '_mstar',
                main_tracer + '_vr_inf', main_tracer + '_r_t', main_tracer + '_f_mol']

fig = multi_panel_plot(x_quantities, y_quantities,
                       parameters, limits, labels, different_symbols=parameters[main_tracer + '_r_bar'],
                       background_cloud=['muse_mass', 'alma'], edge_colours=edge_colours,
                       c_quantities=main_tracer + '_incl',
                       cax_loc=[0.85, 0.1, 0.0125, 0.8], use_errors=True)

# fig = multi_panel_plot(x_quantities, y_quantities,
#                        parameters, limits, labels,
#                        edge_colours=edge_colours, c_quantities=None)

# plt.show()

plt.savefig(plot_folder + 'correlations_global_parameters_' + main_tracer + '.png',
            bbox_inches='tight')
plt.savefig(plot_folder + 'correlations_global_parameters_' + main_tracer + '.pdf',
            bbox_inches='tight')

plt.close()

# FINAL PLOT: Corner plot of the previous y-axes against each other

fig = corner(y_quantities, parameters, limits, labels, different_symbols=parameters[main_tracer + '_r_bar'],
             use_errors=True)

plt.savefig(plot_folder + 'correlations_corner_' + main_tracer + '.png',
            bbox_inches='tight')
plt.savefig(plot_folder + 'correlations_corner_' + main_tracer + '.pdf',
            bbox_inches='tight')

# plt.show()

plt.close()

print('Complete! Took %.2fs' % (time.time() - start))
