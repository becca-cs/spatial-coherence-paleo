# spatial-coherence-paleo

This page is in progress.

Accompanying code and data for "Statistical fingerprints of forced and unforced variability reveal inconsistencies between marine proxies and climate models on multi-decadal to millennial timescales" (Cleveland Stout et al, submitted)

Included here is a zip file "proxies" with proxy data, as well as a "tools" package.

Tools contains:

1. pmtmLS_package, a package for computing the Lomb-Scargle multitaper.
- multitaper_LS_package : contains functions to generate the MTLS; the main function is pmtmLS
- generate_timeseries.m : generates timeseries with a given power spectrum
- testMTLS.m : script for calculating MT and MTLS on a synthetic power spectrum, generated using generate_timeseries

2. sedproxy, a marine sediment proxy-system model package
- age_depth_model.m : age-depth sediment model
- plot_age_depth_example.m : example implementation of age_depth_model.m
- test_ClimToProxyClim.m : example of Sedproxy implementation (see https://github.com/EarthSystemDiagnostics/sedproxy as well for implementation in R)
