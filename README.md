# spatial-coherence-paleo

Accompanying code and data for "Statistical fingerprints of forced and unforced variability reveal inconsistencies between marine proxies and climate models on multi-decadal to millennial timescales" (Cleveland Stout et al, in revision at Paleoceanography and Paleoclimatology)

Included here are two proxy datasets ("MgCa.mat" and "uk37.mat"), as well as a "tools" package.

Tools contains:

1. pmtmLS_package, a package for computing the Lomb-Scargle multitaper.
- multitaper_LS_package : contains functions to generate the MTLS; the main function is pmtmLS
- generate_timeseries.m : generates timeseries with a given power spectrum
- testMTLS.m : script for calculating MT and MTLS on a synthetic power spectrum, generated using generate_timeseries

2. age_depth_model, which produces synthetic age-depth profiles based on BACON (Blauuw and Christen, 2011)
- age_depth_model.m : age-depth sediment model
- plot_age_depth_example.m : example implementation of age_depth_model.m

3. sedproxy, a marine sediment proxy-system model package. This package is modified from the original R
- test_ClimToProxyClim.m : example of Sedproxy implementation (see https://github.com/EarthSystemDiagnostics/sedproxy as well for implementation in R; Dolman et al., 2018)
