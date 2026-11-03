#!/usr/bin/env python3
"""Check the temperature of the intergalactic medium evolved by `intergalacticMediumStateEvolve` before reionization.

Before the ionizing background switches on, the temperature of the intergalactic medium is set by adiabatic expansion,
Compton coupling to the cosmic microwave background, and the sharing of energy with electrons as their number changes:

    dT/dt = -2 H T + (8 σ_T a_r T_γ⁴ / 3 m_e c) (n_e/n_tot) (T_γ - T) - T (dN/dt)/N,

where N = n_tot a³ (only the number of electrons per comoving volume changes, so dN/dt / N = (dn_e/dt + 3 H n_e)/n_tot).
This script runs the IGM evolution model used by `validate-baryonicSuppression.py`, then integrates the equation above
independently, starting from the model's own temperature at its earliest output (z=150) and using the model's own ionization history, and
requires the two temperatures to agree for 8 < z < 140. The atomic heating and cooling terms are negligible over this
range. (Before October 2026 a mis-formed particle-number term, partly offset by spurious recombination cooling of the
neutral gas, made the model's temperature too high by up to 70% at z=8.)

Andrew Benson, Claude
"""
import os
import subprocess
import sys
import xml.etree.ElementTree as ET

import h5py
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

# Relative tolerance on the temperature.
TOLERANCE                   = 1.0e-3
# Range of redshifts over which the temperature is compared.
REDSHIFT_MINIMUM            =   8.0
REDSHIFT_MAXIMUM            = 140.0
# Physical constants (SI). The megaparsec matches Galacticus' definition.
THOMSON_CROSS_SECTION       = 6.6524587321e-29
RADIATION_CONSTANT          = 7.565733250e-16
ELECTRON_MASS               = 9.1093837015e-31
SPEED_LIGHT                 = 2.99792458e8
MEGAPARSEC                  = 3.08567758135e22

# This script runs with the working directory set to testSuite/, while the model is run from the repository root.
os.makedirs("outputs", exist_ok=True)
outputFile = "outputs/igmThermalEvolution.hdf5"
tree = ET.parse("parameters/validate_baryonicSuppression_IGM_evolution.xml")
tree.getroot().find("outputFileName").set("value", "testSuite/"+outputFile)
tree.write("outputs/igmThermalEvolution.xml")
if os.path.exists(outputFile):
    os.remove(outputFile)
status = subprocess.run("cd ..; ./Galacticus.exe testSuite/outputs/igmThermalEvolution.xml", shell=True, capture_output=True)
log = status.stdout.decode() + status.stderr.decode()
if status.returncode != 0 or not os.path.exists(outputFile):
    print("FAILED: the IGM evolution model did not run")
    print(log[-2000:])
    sys.exit(0)
print("SUCCESS: model run")

# Read the IGM state and the cosmological parameters.
with h5py.File(outputFile, "r") as f:
    try:
        igm        = {key: f["igmProperties"][key][:] for key in f["igmProperties"]}
        cosmology  = f["Parameters/cosmologyParameters"].attrs
        H0         = cosmology["HubbleConstant" ]*1.0e3/MEGAPARSEC
        OmegaM     = cosmology["OmegaMatter"    ]
        OmegaL     = cosmology["OmegaDarkEnergy"]
        TCMB0      = cosmology["temperatureCMB" ]
    except KeyError as e:
        print(f"FAILED: expected data missing from the output: {e}")
        sys.exit(0)
redshift           = igm["redshift"]
densityElectron    = igm["densityHydrogen2"]+igm["densityHelium2"]+2.0*igm["densityHelium3"]
densityTotal       = igm["densityHydrogen1"]+igm["densityHydrogen2"]+igm["densityHelium1"]+igm["densityHelium2"]+igm["densityHelium3"]+densityElectron
if redshift.max() <= REDSHIFT_MAXIMUM or redshift.min() >= REDSHIFT_MINIMUM:
    print(f"FAILED: the model's redshifts ({redshift.min():.1f} to {redshift.max():.1f}) do not span the comparison range")
    sys.exit(0)

# Interpolate the ionization history in log(1+z). Comoving electron and total densities are used for the particle number term.
order                  = np.argsort(redshift)
logOnePlusRedshift     = np.log(1.0+redshift[order])
fractionElectron       = interp1d(logOnePlusRedshift,(densityElectron/densityTotal)[order])
logDensityElectronCmvg = interp1d(logOnePlusRedshift,np.log(densityElectron[order]/(1.0+redshift[order])**3),fill_value="extrapolate")

def hubble(z):
    return H0*np.sqrt(OmegaM*(1.0+z)**3+OmegaL)

def rateOfChange(logExpansionFactor, state):
    # Rate of change of temperature with respect to log(a).
    z           = np.exp(-logExpansionFactor)-1.0
    temperature = state[0]
    temperatureCMB = TCMB0*(1.0+z)
    x           = np.log(1.0+z)
    # d ln(n_e a³)/d ln a = -d ln(n_e a³)/d ln(1+z), by finite difference.
    step        = 1.0e-3
    rateNumber  = -(logDensityElectronCmvg(x+step)-logDensityElectronCmvg(x-step))/2.0/step
    rate        = (
        -2.0*temperature
        +8.0*THOMSON_CROSS_SECTION*RADIATION_CONSTANT*temperatureCMB**4/3.0/ELECTRON_MASS/SPEED_LIGHT*fractionElectron(x)*(temperatureCMB-temperature)/hubble(z)
        -temperature*fractionElectron(x)*rateNumber
    )
    return [rate]

compare    = np.sort(redshift[(redshift > REDSHIFT_MINIMUM) & (redshift < REDSHIFT_MAXIMUM)])[::-1]
indexStart = np.argmax(redshift)
solution   = solve_ivp(rateOfChange,[-np.log(1.0+redshift[indexStart]),-np.log(1.0+compare[-1])],[igm["temperature"][indexStart]],t_eval=-np.log(1.0+compare),rtol=1.0e-9,atol=0.0)
if not solution.success:
    print(f"FAILED: the reference integration failed: {solution.message}")
    sys.exit(0)
temperatureReference = solution.y[0]
temperatureModel     = np.array([igm["temperature"][np.argmin(np.abs(redshift-z))] for z in compare])
errorRelative        = np.abs(temperatureModel/temperatureReference-1.0)
worst                = np.argmax(errorRelative)
if errorRelative[worst] > TOLERANCE:
    print(f"FAILED: IGM temperature differs from the adiabatic+Compton solution by {errorRelative[worst]:.2e} at z={compare[worst]:.2f} (model {temperatureModel[worst]:.4f} K, reference {temperatureReference[worst]:.4f} K)")
else:
    print(f"SUCCESS: IGM temperature matches the adiabatic+Compton solution at {compare.size} redshifts in {REDSHIFT_MINIMUM:.0f} < z < {REDSHIFT_MAXIMUM:.0f} (maximum relative difference {errorRelative[worst]:.2e})")
sys.exit(0)
