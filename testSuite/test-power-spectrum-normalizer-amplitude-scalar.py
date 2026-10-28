#!/usr/bin/env python3
import subprocess
import os
import sys
import h5py
import numpy as np

# Test normalization of the power spectrum to a given amplitude of the primordial curvature power spectrum, A_s
# (powerSpectrumNormalizerAmplitudeScalar). The model normalizes to A_s=2.1e-9 with a spectral index n_s=0.90, using the CLASS
# transfer function. Its sigma_8 and its power spectrum on large scales are compared to those computed by CLASS (v3.4.1, the
# version used by Galacticus) itself, run with the same cosmological parameters (h=0.7, omega_b=0.02205, omega_cdm=0.12495,
# Omega_k=0, T_cmb=2.72548 K, N_ur=3.044, YHe=0.2486), A_s, and n_s. The spectral index is chosen to differ from the default
# of CLASS (0.9649), so that a normalization which ignored the spectral index of the model would fail (by about 2.6% in
# sigma_8).
# Andrew Benson (26-October-2026)
pathParameters = os.path.abspath("parameters/powerSpectrumNormalizerAmplitudeScalar.xml")
pathLog        = os.path.abspath("outputs/powerSpectrumNormalizerAmplitudeScalar.log"   )
pathOutput     = os.path.abspath("outputs/powerSpectrumNormalizerAmplitudeScalar.hdf5"  )

# Reference values computed by CLASS.
sigma8CLASS            = 0.832232     # σ₈.
wavenumberCLASS        = 1.0e-2       # Wavenumber (Mpc⁻¹) at which the power spectrum is compared.
powerSpectrumCLASS     = 8.453276e+04 # P(k) (Mpc³) at that wavenumber, at z=0.
toleranceRelative      = 2.0e-3

# Run the model and check for completion.
os.makedirs("outputs",exist_ok=True)
with open(pathLog,"w") as log:
    status = subprocess.run(f"cd ..; ./Galacticus.exe {pathParameters}",stdout=log,stderr=log,shell=True)
if status.returncode != 0:
    print("FAILED: amplitudeScalar power spectrum normalizer model run:")
    subprocess.run(f"cat {pathLog}",shell=True)
    sys.exit(0)
status = subprocess.run(f"grep -q -i -e fatal -e aborted -e \"Galacticus experienced an error in the GSL library\" {pathLog}",shell=True)
if status.returncode == 0:
    print("FAILED: amplitudeScalar power spectrum normalizer model run (errors):")
    subprocess.run(f"cat {pathLog}",shell=True)
    sys.exit(0)
print("SUCCESS: amplitudeScalar power spectrum normalizer model run")

# Extract the power spectrum and σ(M).
model         = h5py.File(pathOutput,'r')
output        = model['Outputs/Output1']
wavenumber    = output['wavenumber'   ][:]
powerSpectrum = output['powerSpectrum'][:]
mass          = output['mass'         ][:]
sigma         = output['sigma'        ][:]
model.close()

# Find σ₈ as σ(M) at the mass contained within a sphere of radius 8/h Mpc.
HubbleConstant     = 0.70
OmegaMatter        = 0.30
densityCritical    = 2.77536627e11*HubbleConstant**2 # M☉ Mpc⁻³.
mass8              = 4.0/3.0*np.pi*(8.0/HubbleConstant)**3*OmegaMatter*densityCritical
order              = np.argsort(mass)
sigma8             = np.exp(np.interp(np.log(mass8          ),np.log(mass      [order]),np.log(sigma        [order])))
powerSpectrumModel = np.exp(np.interp(np.log(wavenumberCLASS),np.log(wavenumber       ),np.log(powerSpectrum       )))

for name, value, reference in (
        ("σ₈"                                , sigma8            , sigma8CLASS       ),
        (f"P(k={wavenumberCLASS} Mpc⁻¹, z=0)", powerSpectrumModel, powerSpectrumCLASS)
):
    difference = value/reference-1.0
    if abs(difference) < toleranceRelative:
        print(f"SUCCESS: {name} = {value:.6e} agrees with CLASS ({reference:.6e}; relative difference {difference:+.2e})")
    else:
        print(f"FAILED: {name} = {value:.6e} disagrees with CLASS ({reference:.6e}; relative difference {difference:+.2e})")
sys.exit(0)
