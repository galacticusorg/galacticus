#!/usr/bin/env python3
import subprocess
import sys
import numpy as np
import h5py

# Test that the "thin disk" AGN SED model produces SEDs that correctly
# integrate to the bolometric luminosity of the accretion disk.
# Andrew Benson (16-April-2026)

# Make output directory.
subprocess.run("mkdir -p outputs/", shell=True)

# Run the model.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/thinDiskAGNSED.xml", shell=True)
if status.returncode != 0:
    print("FAIL: AGN SED model failed to run")
    sys.exit(0)

# Read the spectra and associated black hole data.
model               = h5py.File('outputs/thinDiskAGNSED.hdf5','r')
nodes               = model['Outputs/Output1/nodeData']
massBlackHole       = nodes['blackHoleMass'                ][:]
rateBlackHole       = nodes['massAccretionRateBlackHoles'  ][:]
efficiencyBlackHole = nodes['radiativeEfficiencyBlackHoles'][:]
wavelength          = nodes['agnSEDColumnValues'           ][:]
sedBlackHole        = nodes['agnSED'                       ][:]

# Define physical constants (SI units).
massSolar             = 1.99e+30
gigaYear              = 3.16e+16
speedOfLight          = 3.00e+08
gravitationalConstant = 6.67e-11
thomsonCrossSection   = 6.65e-29
massHydrogenAtom      = 1.66e-27
luminositySolar       = 3.83e+26
angstrom              = 1.00e-10

# Compute the Eddington ratio.
rateEddington = 4.0*np.pi*gravitationalConstant*massBlackHole*massHydrogenAtom*gigaYear/thomsonCrossSection/speedOfLight
lambdaEdd     = rateBlackHole/rateEddington

# Compute bolometric luminosity.
luminosityBolometric = efficiencyBlackHole*(rateBlackHole*massSolar/gigaYear)*speedOfLight**2/luminositySolar

# Compute frequency grid.
frequency      = speedOfLight/(wavelength*angstrom)
deltaFrequency = -np.diff(frequency)

# Select only galaxies with significant accretion rates.
selection = (lambdaEdd > 0.01) & (rateBlackHole > 1e6)

# Integrate SEDs across frequency.
luminositySED = luminosityBolometric[selection]
for i in range(np.count_nonzero(selection)):
    luminositySED[i] = np.sum(sedBlackHole[selection][i][0:-1]*deltaFrequency)

# Compute the ratio of the SED luminosity to the bolometric luminosity.
ratio = luminositySED/luminosityBolometric[selection]

# Test that all ratios are close to 1.
if np.allclose(ratio,1.0,rtol=0.0,atol=0.03):
    print("SUCCESS: AGN SEDs integrate to bolometric luminosity")
else:
    print("FAIL: AGN SEDs do not integrate to bolometric luminosity")
