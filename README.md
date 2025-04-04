# TPM wrapper
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Developer
- Jin BENIYAMA [mail](mailto:jinbeniyama@gmail.com)


## Overview
Here are useful scripts to do TPMs as well as to plot results of TPMs.
Under heavy development by J.B.


### Examples
```
# Collect observations in a single file (obs.dat) as below (by hand!)
``
jd wavelength flux fluxerr code cflag memo
2452545.9060276826 8.12 6.290097664831847 0.07037871513098569 568 999 UKIRT2002
2452545.9060276826 8.369 6.98550199651787 0.04672576586299578 568 999 UKIRT2002
``

# Or you can make obs.dat without observation.
make_obsdat.py 2460457.7236111113 2460457.8194444445 0.0007 --wavelength 8.7 10.7 11.7 --code 309 --out obs20240527.dat


# Make obs and ephem file from observations of (433) Eros with light-time correction
make_obseph.py 433 (obs.dat) --out_obs obs.txt --out_eph eph.txt --ltcor
# Make obs and ephem file from observations of (433) Eros without light-time correction
make_obseph.py 433 (obs.dat) --out_obs obs.txt --out_eph eph.txt

# Make obs and ephem file for simulation
# P=60 s, t_sample=1 s, 2 rotation phase, alpha = 60 deg
make_obseph_sim.py --obs obs_alpha60_20240515.txt --eph eph_alpha60_20240515.txt --rotP_s 60.0 --t_sample 1 --N_rot 2 --alpha 60

# Plot results of TPM at once (TI vs. reduced chi square)
plot_tpm_brute.py tpmout* -x TI --reduce

# Plot results of TPM at once (A vs. reduced chi square)
plot_tpm_brute.py tpmout* -x A --reduce

# Plot thermal lightcurve
plot_thermal_flux.py (output of tpm)

# Blend thermal fluxes
blend_tpm_result.py tpmout*

# Plot 3-d chi-square map
plot_blended_result.py results.txt
```

## Dependencies
This library is depending on `NumPy`, `SciPy`, `SEP`, `Astropy` 
and `Astroquery`.
Scripts are developed on `Python 3.7.10`, `NumPy 1.19.2`, `SciPy 1.6.1`,
`SEP 1.0.3`, `Astropy 4.2` and `Astroquery 0.4.1`.


## LICENCE
This software is released under the MIT License.
