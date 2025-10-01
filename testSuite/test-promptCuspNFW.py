#!/usr/bin/env python3
# Test that calculations of prompt cusp-NFW profile properties match the results from Sten Delos' implementation.
import sys
sys.path.append("../cusp-halo-relation")
import numpy as np
import h5py
import subprocess
import os
from cusp_halo_relation import cuspNFW

# Create output path.
try:
    os.mkdir("outputs")
except FileExistsError:
    pass

# Run the models.
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/parameters/testPromptCuspNFW.xml",shell=True)
if status.returncode != 0:
   print("FAILED: model failed to run"  )
   sys.exit()

# Read require data from the model.
model                 = h5py.File('outputs/testPromptCuspNFW.hdf5','r')
nodes                 = model['Outputs/Output1/nodeData']
massVirial            = nodes['basicMass'                                 ][:]
radiusVirial          = nodes['darkMatterOnlyRadiusVirialLastDefined'     ][:]
radiusMinus2          = nodes['darkMatterProfileScale'                    ][:]
amplitudeCusp         = nodes['darkMatterProfilePromptCuspAmplitude'      ][:]
radiusScale           = nodes['darkMatterProfilePromptCuspNFWRadiusScale' ][:]
densityScale          = nodes['darkMatterProfilePromptCuspNFWDensityScale'][:]

# Compute derived quantities.
densityVirial         = massVirial/(4*np.pi/3*radiusVirial**3)
concentration         = radiusVirial/radiusMinus2

# Compute reference quantities.
radiusScaleReference  = np.zeros_like(amplitudeCusp)
densityScaleReference = np.zeros_like(amplitudeCusp)
for i in range(len(amplitudeCusp)):
  radiusScaleReference[i],densityScaleReference[i] = cuspNFW.scale_from_c(concentration[i],massVirial[i],amplitudeCusp[i],densityVirial[i],cmin_error=False)

# Evaluate maximum errors in scale radius and density.
errorFractionalRadiusScale       = np.abs( radiusScale- radiusScaleReference)/ radiusScaleReference
errorFractionalDensityScale      = np.abs(densityScale-densityScaleReference)/densityScaleReference
errorFractionalRadiusScaleWorst  = np.max(errorFractionalRadiusScale )
errorFractionalDensityScaleWorst = np.max(errorFractionalDensityScale)

# Validate the results.
if errorFractionalRadiusScaleWorst > 1.0e-3 or errorFractionalDensityScaleWorst > 1.0e-3:
    print("FAILED: results do not agree with reference calculation")
else:
    print("SUCCESS: results do agree with reference calculation")
