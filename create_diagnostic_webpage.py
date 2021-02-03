# -*- coding: utf-8 -*-
"""
Create the webpage for diagnostics

@author: Tom Williams
"""

import os

from make_html_classpage import make_classification_page
from vars import phangs_folder, alma_version, alma_plot, muse_plot, muse_version, alma_galaxies

os.chdir(phangs_folder)

webpage_dir = '/Users/williams/PycharmProjects/pattern_speed_webpage'

emission_types = ['alma', 'mass', 'ha']

for emission_type in emission_types:

    html_folder = os.path.join(webpage_dir, 'images')
    if not os.path.exists(html_folder):
        os.makedirs(html_folder)

    # Copy the files we need (just the <v>/<x> plots)

    for galaxy in alma_galaxies:

        galaxy = galaxy.upper()

        if 'alma' in emission_type:
            file_in_name = os.path.join(alma_plot, alma_version, 'bmask', galaxy + '_bmask_alma.png')
        else:
            file_in_name = os.path.join(muse_plot, muse_version, emission_type + '_smask_bmask',
                                        galaxy + '_' + emission_type + '_smask_bmask_muse.png')

        file_out_name = os.path.join(html_folder, galaxy + '_' + emission_type + '.png')

        if not os.path.exists(file_out_name) and os.path.exists(file_in_name):
            os.system('cp ' + file_in_name + ' ' + file_out_name)

        # Now, create the webpage

        make_classification_page(folder=webpage_dir,
                                 fig_folder='images',
                                 prefix="", suffix="",
                                 nbatches=2, person_name='Tom Williams', email='williams@mpia.de',
                                 list_classes=["Single-well-definedP", "ClearMultiplebutClean", "PoorFit",
                                               "InsuffQuality"],
                                 list_explanations=[
                                     "Single, well defined pattern speeds: no issues (integrals have converged, "
                                     "reasonably stable with slit width, sufficiently high covering factor",
                                     'Clear multiple pattern speeds visible in the <v><x> plot, but otherwise the fit '
                                     'would be a quality flag (1)',
                                     "Poor fit: integral hasn't converged, or there is some issue in the data that "
                                     "causes us to distrust this pattern speed",
                                     'Data of insufficient quality to calculate a reliable pattern speed. In this '
                                     'case, a more reliable pattern speed may be possible with higher resolution or '
                                     'deeper data. '
                                 ],
                                 title='Classification of Pattern Speed for ALMA/MUSE',
                                 subtitle='Classification of Pattern Speed for ALMA/MUSE'
                                 )
