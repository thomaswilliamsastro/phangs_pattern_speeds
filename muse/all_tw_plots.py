# -*- coding: utf-8 -*-
"""
Create a tex file that'll plot all the TW MUSE stellar mass plots

@author: Tom Williams
"""

import os
import shutil

from astropy.table import Table

from vars import phangs_folder, muse_version, muse_galaxies, muse_plot, pattern_speed_version, output_folder

os.chdir(phangs_folder)

pattern_speed_table = Table.read(output_folder + 'pattern_speed_table_' + pattern_speed_version + '.fits')

f = open('pattern_speeds_output/all_tw_plots.tex', 'w+')
f_arxiv = open('pattern_speeds_output/all_tw_plots_arxiv.tex', 'w+')

if not os.path.exists('pattern_speeds_output/all_muse_tw'):
    os.makedirs('pattern_speeds_output/all_muse_tw')

f.write('\\figsetstart\n')
f.write('\\figsetnum{1}\n')
f.write('\\figsetgrpstart\n')

for i, galaxy in enumerate(muse_galaxies):

    table_row = pattern_speed_table[pattern_speed_table['GALAXY'] == galaxy.upper()]
    q = int(table_row['OM_P_MUSE_MASS_QUAL'][0])

    if i == 0:
        example_q = q

    muse_file = os.path.join(muse_plot, muse_version, 'mass_smask_bmask', galaxy + '_mass_smask_bmask_muse.pdf')

    shutil.copy(muse_file, 'pattern_speeds_output/all_muse_tw/' + galaxy.upper() + '_mass.pdf')

    f_arxiv.write('\\begin{figure*}[t]\n')
    f_arxiv.write('\plotone{%s_mass.pdf}\n' % galaxy.upper())
    f_arxiv.write('\caption{As Fig. \\ref{fig:ngc3351_tw_integral}, but for %s. For this galaxy , $Q=%s$. \label{fig:app_%s}}\n'
            % (galaxy.upper(), q, galaxy))
    f_arxiv.write('\end{figure*}\n')
    f_arxiv.write('\n')

    f.write('\\figsetgrpnum{B.%d}\n' % (i+1))
    f.write('\\figsetgrptitle{%s}\n' % galaxy.upper())
    f.write('\\figsetplot{%s_mass.pdf}\n' % galaxy.upper())
    f.write('\\figsetgrpnote{As Fig. \\ref{fig:ngc3351_tw_integral}, but for %s. For this galaxy , $Q=%s$.}\n'
            % (galaxy.upper(), q))
    f.write('\\figsetgrpend\n')

f_arxiv.close()

f.write('\\figsetgrpend')

# Finally, include an example Fig.

f.write('\\begin{figure*}[t]\n')
f.write('\plotone{NGC3351_mass_smask_bmask_muse.pdf}\n')
f.write('\caption{{\it Top left}: Stellar mass surface density map of NGC~3351 shown in greyscale, with '
        'Tremaine-Weinberg integral slits of 1\\arcsec~width, oriented parallel to the major axis overlaid. Only one in '
        'every four slits is shown, due to the slit density, and are coloured according to their position along the '
        'kinematic minor axis. For this galaxy, the quality flag, $Q = 1$ (see Sect. \\ref{sec:quality_flagging}). '
        '{\it Top right}: Stellar velocity map for the same galaxy, with a dashed black line showing the kinematic '
        'major axis, passing through the galaxy centre. {\it Bottom}: intensity-weighted velocity ($\langle v \\rangle$) '
        'versus intensity-weighted position ($\langle x \\rangle$) for each of the slits (the colour corresponds to the '
        'slit colour in the above top-left panel). The black line shows the best fit, and the grey shaded region the '
        'errors on the fit (in this case, this region is extremely small). One point has an extremely large uncertainty '
        'in this panel, and the error bar extends across the entire range of $\langle v \\rangle$ shown. '
        'The complete figure set (19 images) is available in the online journal. \label{fig:ngc3351_tw_integral}}')
f.write('\end{figure*}\n')
f.write('\n')

f.close()

print('Complete!')
