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
  sys.exit(0)

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
countHalos            = len(amplitudeCusp)

# We validate the scale parameters by evaluating the residuals of the two equations which define them, rather than by comparing
# rₛ and ρₛ directly against a reference solve. A direct comparison of rₛ is ill-conditioned: the mass equation becomes
# insensitive to rₛ as a halo approaches the minimum-concentration boundary (∂M/∂rₛ → 0 there), so a solution of fixed quality
# shows an error in rₛ which grows without bound as that boundary is approached. Halos sitting a fractional margin of ~2×10⁻⁴
# from the boundary amplify the residual by a factor of ~150. The residuals used below are, by contrast, uniform across the
# full range of margins.
radiusVirialReference = np.array([cuspNFW.R_from_M(massVirial[i],densityVirial[i]) for i in range(countHalos)])

# Residual of the mass equation: the mass of the cusp-NFW profile enclosed within the virial radius must equal the halo mass.
# This is the defining constraint on (rₛ,ρₛ).
massEnclosed          = np.array([cuspNFW.mass(radiusVirialReference[i],radiusScale[i],densityScale[i],amplitudeCusp[i]) for i in range(countHalos)])
residualMass          = np.abs(np.log(massEnclosed/massVirial))

# Residual of the geometric relation between the scale radius and r₋₂. Halos which have been pushed to the minimum-concentration
# limit are excluded here: for those Galacticus sets rₛ from the limiting value of the y-parameter instead of from this
# relation, so the relation does not apply to them. (Their mass-equation residuals are tested above along with all others.)
concentrationMinimum  = np.array([cuspNFW.c_min(massVirial[i],amplitudeCusp[i],densityVirial[i]) for i in range(countHalos)])
atConcentrationLimit  = concentration < concentrationMinimum
radiusScaleFromMinus2 = np.array([cuspNFW.rs_from_r2(radiusVirialReference[i]/concentration[i],densityScale[i],amplitudeCusp[i]) for i in range(countHalos)])
residualGeometry      = np.abs(radiusScaleFromMinus2-radiusScale)/radiusScale

residualMassWorst     = np.max(residualMass                          )
residualGeometryWorst = np.max(residualGeometry[~atConcentrationLimit])

# Validate the results.
if residualMassWorst > 2.0e-4 or residualGeometryWorst > 1.0e-10:
    print("FAILED: results do not agree with reference calculation (worst mass residual = {:.3e}, worst geometry residual = {:.3e})".format(residualMassWorst,residualGeometryWorst))
else:
    print("SUCCESS: results do agree with reference calculation")
