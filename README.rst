#####################
PHANGS Pattern Speeds
#####################

.. note::
   TODO

   * Categorise pattern speeds into 1) looks good; 2) looks ok; 3) don't trust

The code in this repository calculates the bar pattern speed for the PHANGS MUSE and PHANGS ALMA galaxies using the
Tremaine-Weinberg method. This document provides an overview of what each file does. The repository is structured so
that any code specific to MUSE is contained within the MUSE folder, and the same with ALMA. Any code that is used by
both is in this top folder.

---------------------------
``alma_muse_tw_comparison``
---------------------------

Produces a plot comparing the measured pattern speeds between ALMA and MUSE. Well-constrained pattern speeds are filled
circles, and the bad fits are unfilled.

--------------
``bootstraps``
--------------

Contains the functions that allow for bootstrapping errors in centering and position angle for calculating pattern
speeds.

-------------------------
``create_diagnostic_doc``
-------------------------

Produces a PDF file for diagnostic purposes. Each page will contain one galaxy, showing the TW integrals for emission
masked beyond the bar and not, as well as the slit lengths and widths effect on the recovered pattern speeds.

-----------------------------
``literature_pattern_speeds``
-----------------------------

Produces a plot comparing the measured pattern speeds to published literature values (homogenised to our distance and
inclination values). Well-constrained pattern speeds are filled circles, and the bad fits are unfilled.

--------------
``make_table``
--------------

Produces a .fits and .tex table of pattern speeds and corotation radii.

----------------
``ps_functions``
----------------

Contains useful functions for getting at the PHANGS master table and pre-processing the data for the TW integrals.

---------------------------
``quality_flag_morph_type``
---------------------------

Produces a plot showing the distribution of well-measured pattern speeds (and the underlying sample) against
morphological type.

-----------
``r_corot``
-----------

Using the measure pattern speeds, calculates an arbitrary number of corotation radii using the fitting rotation curves
from Lang+ subm.

-----------------------------
``radial_profile_comparison``
-----------------------------

Produces a plot of corotation radii versus radially averaged quantities (measured at 200pc by Jiayi Sun).

--------------------------------
``reduction_version_comparison``
--------------------------------

Produces a plot comparing the measured pattern speeds for different ALMA/MUSE data reduction versions.

--------
``vars``
--------

Contains paths to a number of folders, and a number of frequently used variables in this work.

====
ALMA
====

-----------------------
``bar_mask_comparison``
-----------------------

Produces a plot comparing the measured pattern speeds if emission beyond the bar is masked, or not.

------------
``get_data``
------------

Connects to the PHANGS server at MPIA and downloads the current latest reductions of the ALMA mom0, mom1 and associated
error maps at the highest resolution (priority is '12m+7m+tp', '12m+7m', '7m+tp', then finally '7m'). This will work if
you're at MPIA, but may need some editing if you're downloading from a different server.

------------------
``plot_tw_speeds``
------------------

Produces a plot of the different pattern speeds for each galaxy. As there are a lot of ALMA sources, also produces a KDE
plot for the overall distribution of pattern speeds.

-----------------
``sl_comparison``
-----------------

Produces a diagnostic plot for the ALMA data, comparing the recovered pattern speed for a variety of different slit
radii. The recovered pattern speed tend to level out at large radii, so to aid convergence I allow the slit to be as
long as possible.

-----------------
``sw_comparison``
-----------------

Produces a diagnostic plot for the ALMA data, comparing the recovered pattern speed for a variety of different slit
widths. The recovered pattern speed tends to be constant up to around 7 or so arcsec, so to make the most of the
resolution of the data I leave it at 1 arcsec.

-----------
``tw_alma``
-----------

Sets up ALMA maps to be fed into ``bootstraps`` to calculate a pattern speed. Contains a number of options which have
been used in testing.

====
MUSE
====

-----------------------
``coverage_comparison``
-----------------------

Produces a diagnostic plot for the MUSE data, comparing the recovered pattern speed depending on the covering factor of
the H-alpha maps (based on a per-pixel S/N cutoff). Given the high S/N of the data across the disc, this doesn't affect
the recovered pattern speed.

-----------------------
``emission_comparison``
-----------------------

Produces a diagnostic plot for the MUSE data, comparing the recovered pattern speed for a variety of different
flux/velocity tracers (stellar mass, white light, H-alpha), as well as the effect of masking stars and regions beyond
the bar.

------------
``get_data``
------------

Connects to the PHANGS server at MPIA and downloads the current latest reductions of the output MAPS files (the
multi-extension .fits files containing all the various parameter maps). This will work if you're at MPIA, but may need
some editing if you're downloading from a different server.

------------------
``plot_tw_speeds``
------------------

Produces a comparison plot of the pattern speeds for all the MUSE galaxies.

-----------------
``sl_comparison``
-----------------

Produces a diagnostic plot for the MUSE data, comparing the recovered pattern speed for a variety of different slit
radii. The recovered pattern speed tend to level out at large radii, so to aid convergence I allow the slit to be as
long as possible.

-----------------
``sw_comparison``
-----------------

Produces a diagnostic plot for the MUSE data, comparing the recovered pattern speed for a variety of different slit
widths. The recovered pattern speed tends to be constant up to around 7 or so arcsec, so to make the most of the
resolution of the data I leave it at 1 arcsec.

-----------
``tw_muse``
-----------

Sets up MUSE maps to be fed into ``bootstraps`` to calculate a pattern speed. Contains a number of options which have
been used in testing.

=======
Changes
=======

----------
2020/03/04
----------

* Code uses data reduction versions throughout.
* Comparison between different data reductions is now possible.
* Attempt at fitting multiple pattern speeds added -- seems like this doesn't work :(
* Corotation radii are now calculated.
* H-alpha pattern speeds are also calculated.
* Diagnostic documents produced, which allow for (hopefully) easy spotting of bad pattern speeds.
* Initial pass of quality flagging done (TW), and code updated to reflect the bad fits.
* v0.1 of the table is now available, containing the pattern speeds, corotation radii and initial quality flags.
* Comparison to published literature values included.
* Initial results using radial profiles added.
* Many other smaller things.

----------
2019/12/13
----------

* Modified script to produce overview plot for ALMA. Since there are a lot of galaxies, also a KDE plot
  (``alma/plot_tw_speeds``)
* Modified scripts to test slit widths/lengths for the ALMA test case (``alma/sw_comparison`` and
  ``alma/sl_comparison``)
* Modified MUSE data prep script for ALMA (``alma/tw_alma``)

----------
2019/12/12
----------

* Added in diagnostics for different slit lengths for MUSE (``muse/sl_comparison``)
* Added in diagnostics for different slit widths for MUSE (``muse/sw_comparison``)
* Added in diagnostics for different flux/velocity maps, and various types of mask for the MUSE data
  (``muse/emission_comparison``)
* Added script to bootstrap errors in position angle and centering, to produce reasonable error estimates on pattern
  speed (``bootstraps``)
* Added script to prepare data for bootstrapping, with plenty of testing options for the MUSE data (``muse/tw_muse``)