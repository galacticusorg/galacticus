# SIDM parametric model test data

Reference data used by `source/tests/SIDM_parametric_model.F90`, extracted from an N-body simulation.

- `data_799_cdm_NFW.txt` — CDM NFW halo mass accretion history. One header line, then rows of:
  `expansionFactor  virialMass[Msun/h]  virialRadius[kpc/h,comoving]  radiusVmax[kpc/h,comoving]`
- `data_799_cdm.txt`     — CDM halo properties. One header line, then columns of which column 7 is the
  maximum circular velocity Vmax [km/s].

Both files have the same number of data rows, tabulated from earliest to latest time.
