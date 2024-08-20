# spatial-coherence-paleo
Accompanying code and data for "Statistical fingerprints of forced and unforced variability reveal inconsistencies between marine proxies and climate models on multi-decadal to millennial timescales" (Cleveland Stout et al, submitted)

Included here is a zip file "proxies" with proxy data, as well as a "tools" package.

Tools contains:
- pmtmLS_package : a package for computing the Lomb-Scargle multitaper.
-     multitaper_LS_package : contains functions to generate the MTLS; the main function is pmtmLS
-     generate_timeseries.m : generates timeseries with a given power spectrum
-     testMTLS.m : script for calculating MT and MTLS on a synthetic power spectrum, generated using generate_timeseries
- 
