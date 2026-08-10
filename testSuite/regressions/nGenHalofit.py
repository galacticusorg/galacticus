#!/usr/bin/env python3
import h5py
import numpy as np
import subprocess
import sys

# Run a power spectrum task using the NGenHalofit nonlinear power spectrum twice. This tests that a
# pre-computed tabulation of the ratio of nonlinear to linear power is correctly re-read on subsequent
# runs. Then check that the resulting nonlinear power spectrum behaves as expected relative to the
# linear power spectrum.
# Andrew Benson (08-August-2026)

# Run the model and check for successful completion.
for i in range(2):
    status = subprocess.run("./Galacticus.exe testSuite/regressions/nGenHalofit.xml", shell=True)
    if status.returncode != 0:
        print("FAILED: NGenHalofit model failed to complete")
        sys.exit(0)

# Read the resulting power spectra. Output redshifts are computed from the corresponding output times, so
# will not exactly equal the requested values - match them to within a tolerance instead.
model      = h5py.File("testSuite/outputs/regressions/powerSpectrumNGenHalofit.hdf5","r")
ratio      = {}
wavenumber = None
for name in model["Outputs"]:
    group = model["Outputs"][name]
    if "powerSpectrumNonlinear" not in group:
        print("FAILED: NGenHalofit model did not output a nonlinear power spectrum")
        sys.exit(0)
    powerLinear    = group["powerSpectrum"         ][:]
    powerNonlinear = group["powerSpectrumNonlinear"][:]
    wavenumber     = group["wavenumber"            ][:]
    for redshift in (0.0, 1.0):
        if np.abs(group.attrs["outputRedshift"]-redshift) < 1.0e-6:
            ratio[redshift] = powerNonlinear/powerLinear

for redshift in (0.0, 1.0):
    if redshift not in ratio:
        print("FAILED: NGenHalofit model did not output at z="+str(redshift))
        sys.exit(0)

# The ratio of nonlinear to linear power must be close to unity on large scales, and must boost power
# on small scales. Wavenumbers here are in units of Mpc¯¹; NGenHalofit is calibrated against
# simulations over 0.05 h/Mpc ≤ k ≤ 10 h/Mpc. The "small" range is kept inside the tabulated range of
# wavenumbers (which extends to 100 h/Mpc) so that no extrapolation is tested here.
large = wavenumber < 1.0e-2
small = (wavenumber > 1.0e+0) & (wavenumber < 5.0e+1)
if not np.any(large) or not np.any(small):
    print("FAILED: NGenHalofit power spectrum does not span the wavenumbers to be tested")
    sys.exit(0)
for redshift in (0.0, 1.0):
    if np.any(np.abs(ratio[redshift][large]-1.0) > 0.05):
        print("FAILED: NGenHalofit nonlinear power differs from linear power on large scales at z="+str(redshift))
        sys.exit(0)
    if np.any(ratio[redshift][small] < 1.0):
        print("FAILED: NGenHalofit does not boost power on small scales at z="+str(redshift))
        sys.exit(0)

# Nonlinear evolution must be less advanced at higher redshift, so the ratio of nonlinear to linear
# power on small scales must be smaller at z=1 than at z=0.
if np.any(ratio[1.0][small] >= ratio[0.0][small]):
    print("FAILED: NGenHalofit nonlinear enhancement does not decrease with increasing redshift")
    sys.exit(0)

print("SUCCESS: NGenHalofit")
