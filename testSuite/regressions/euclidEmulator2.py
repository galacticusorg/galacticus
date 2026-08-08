#!/usr/bin/env python3
import h5py
import numpy as np
import subprocess
import sys

# Run a power spectrum task using the EuclidEmulator2 nonlinear power spectrum twice. This tests that
# a pre-computed nonlinear correction factor file is correctly re-read on subsequent runs. Then check
# that the resulting nonlinear power spectrum behaves as expected relative to the linear power
# spectrum.
# Andrew Benson (07-August-2026)

# Run the model and check for successful completion.
for i in range(2):
    status = subprocess.run("./Galacticus.exe testSuite/regressions/euclidEmulator2.xml", shell=True)
    if status.returncode != 0:
        print("FAILED: euclidEmulator2 model failed to complete")
        sys.exit(0)

# Read the resulting power spectra.
model  = h5py.File("testSuite/outputs/regressions/powerSpectrumEuclidEmulator2.hdf5","r")
group  = model["Outputs/Output1"]
if "powerSpectrumNonlinear" not in group:
    print("FAILED: euclidEmulator2 model did not output a nonlinear power spectrum")
    sys.exit(0)
wavenumber     = group["wavenumber"            ][:]
powerLinear    = group["powerSpectrum"         ][:]
powerNonlinear = group["powerSpectrumNonlinear"][:]

# The nonlinear correction factor must be close to unity on large scales, and must boost power on
# small scales. Wavenumbers here are in units of Mpc¯¹; EuclidEmulator2 emulates 8.73e-3 h/Mpc ≤ k ≤
# 9.41 h/Mpc.
boost = powerNonlinear/powerLinear
large = wavenumber < 1.0e-2
small = wavenumber > 1.0e+0
if not np.any(large) or not np.any(small):
    print("FAILED: euclidEmulator2 power spectrum does not span the wavenumbers to be tested")
    sys.exit(0)
if np.any(np.abs(boost[large]-1.0) > 0.05):
    print("FAILED: euclidEmulator2 nonlinear correction factor differs from unity on large scales")
    sys.exit(0)
if np.any(boost[small] < 1.0):
    print("FAILED: euclidEmulator2 nonlinear correction factor does not boost power on small scales")
    sys.exit(0)

print("SUCCESS: euclidEmulator2")
