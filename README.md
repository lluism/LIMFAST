# 21cmFAST

An edited version of Rick Mebane's fork of the public, development version of the semi-numerical, cosmological simulation code 21cmFAST (Mesinger & Furlanetto 2007; Mesinger et al. 2011).

# General source models

Functions to use general source models are in Parameter_files/SOURCES.H. To use, set USE_GENERAL_SOURCES = 1 in Parameter_files/ANAL_PARAMS.H (note: SHARP_CUTOFF must be set to 0 as well). For now, the code will always use the functions set in defaultSources(). You can provide your own functions for minMass, fstar, fesc, Nion, fx, and fPopIII which can be functions of both halo mass and redshift. Populations are not explicitly set, so any differences (e.g., Pop II/III) are just set by the properties of these populations in the source functions.
