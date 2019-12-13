#####################
PHANGS Pattern Speeds
#####################

The code in this repository calculates the bar pattern speed for the PHANGS MUSE and PHANGS ALMA galaxies using the
Tremaine-Weinberg method. This document provides an overview of what each file does. The repository is structured so
that any code specific to MUSE is contained within the MUSE folder, and the same with ALMA. Any code that is used by
both is in this top folder.

====
MUSE
====

-----------------------
``emission_comparison``
-----------------------

Produces a diagnostic plot for the MUSE data, comparing the recovered pattern speed for a variety of different
flux/velocity tracers (stellar mass, white light, H-alpha), as well as the effect of masking stars and regions beyond
the bar.

Given the resilience of the recovered pattern speed to these various different methodologies, it seems these aren't
particularly important. I've decided to go with stellar mass and velocity, masking stars and regions beyond the bar,
which seems like the most natural choice.

-----------
``folders``
-----------

Contains variables for the folder names used a lot throughout this MUSE work

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
radii (30 arcsec to 110 arcsec in 5 arcsec steps). The recovered pattern speed levels out at around the bar radius, and
then remains constant, so I choose to leave the slit lengths as long as possible.

-----------------
``sw_comparison``
-----------------

Produces a diagnostic plot for the MUSE data, comparing the recovered pattern speed for a variety of different slit
widths (1 arcsec to 10 arcsec in 0.5 arcsec steps). The recovered pattern speed is the same no matter the slit width,
so I leave it at the default (which will help to show up any interesting features and potential multiple pattern speeds.

-----------
``tw_muse``
-----------

Sets up MUSE maps to be fed into ``bootstraps`` to calculate a pattern speed. Contains a number of options which have
been used in testing. The ones it's set to are "optimal", so should be fine to use as-is. Has options for varying slit
length, slit width, and tracer.

====
ALMA
====

-----------
``folders``
-----------

Contains variables for the folder names used a lot throughout this ALMA work

------------
``get_data``
------------

Connects to the PHANGS server at MPIA and downloads the current latest reductions of the ALMA mom0, mom1 and associated
error maps at the highest resolution (priority is '12m+7m+tp', '12m+7m', '7m+tp', then finally '7m'). This will work if
you're at MPIA, but may need some editing if you're downloading from a different server.

-----------------
``sl_comparison``
-----------------

Produces a diagnostic plot for the MUSE data, comparing the recovered pattern speed for a variety of different slit
radii (10 arcsec to 50 arcsec in 5 arcsec steps). The recovered pattern speed levels out at large radii, so much like
the MUSE I don't do any masking of slit lengths.

-----------------
``sw_comparison``
-----------------

Produces a diagnostic plot for the ALMA data, comparing the recovered pattern speed for a variety of different slit
widths (1 arcsec to 10 arcsec in 0.5 arcsec steps). The recovered pattern speed is the same no matter the slit width,
so I leave it at the default (which will help to show up any interesting features and potential multiple pattern speeds.

-----------
``tw_alma``
-----------

Sets up ALMA maps to be fed into ``bootstraps`` to calculate a pattern speed. Contains a number of options which have
been used in testing. The ones it's set to are "optimal", so should be fine to use as-is. Has options for varying slit
length, and slit width.

=======
Changes
=======

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