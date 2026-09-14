#!/usr/bin/env python3
"""Check energy balance in the thermal emission from dust.

Dust must re-emit exactly the luminosity it absorbs. For each of several attenuators, the luminosity absorbed by each
phase of dust is extracted by a `dustAttenuation` extractor (`outputAbsorbed`), and the thermal emission of the dust by
a `SEDDustEmission` extractor heated by the same light. The absorbed luminosity is integrated here, independently of
the emission extractor, and compared with the integral of the emitted spectrum, in total and phase by phase. Heating by
the cosmic microwave background is switched off, since with it the dust also re-emits energy absorbed from the CMB.

A second output, at redshift one, checks that the observed-frame emission on a grid twice as long in wavelength equals
the rest-frame emission bin by bin.

Andrew Benson, with assistance from Claude.
"""
import os
import subprocess
import sys

import h5py
import numpy as np

# Tolerance on energy balance. The absorbed luminosity is integrated by the trapezoidal rule on the grid of the heating
# spectrum, and the emission by summing over resolution elements, so the two agree only to the accuracy of those
# quadratures.
TOLERANCE_ENERGY = 2.0e-3
# Tolerance on the observed-frame emission matching the rest-frame emission, limited by the precision with which the
# expansion factor at the output time reproduces exactly one half.
TOLERANCE_FRAME  = 1.0e-6
# Speed of light in Å s⁻¹, and the Solar luminosity in W as used by Galacticus.
SPEED_OF_LIGHT   = 2.99792458e18
LUMINOSITY_SOLAR = 3.845e26
# Resolution of the emission spectra in the parameter file, and the factor by which the extremes of a resolution
# element differ from its center.
RESOLUTION       = 50.0
FACTOR           = (1.0 + np.sqrt(1.0 + 4.0 * RESOLUTION**2)) / 2.0 / RESOLUTION

ATTENUATORS = {
    "charlotFall2000"           : ("birthCloud", "screenSurfaceDensityMetals"),
    "screenSurfaceDensityMetals": ("screenSurfaceDensityMetals",),
    "atlasFerrara2000"          : ("atlasFerrara2000",),
}

# This script runs with the working directory set to testSuite/, while the model is run from the repository root.
parameterFile = "testSuite/parameters/dustEmissionEnergyBalance.xml"
outputFile    = "testSuite/outputs/dustEmissionEnergyBalance.hdf5"
outputPath    = os.path.join("..", outputFile)

os.makedirs("outputs", exist_ok=True)
if os.path.exists(outputPath):
    os.remove(outputPath)

status = subprocess.run(f"cd ..; ./Galacticus.exe {parameterFile}", shell=True, capture_output=True)
log = status.stdout.decode() + status.stderr.decode()
if status.returncode != 0 or not os.path.exists(outputPath):
    print("FAILED: the dust emission energy balance model did not run")
    print(log[-2000:])
    sys.exit(0)


def columnValues(nodes, name, fallback):
    """Return the wavelengths of the columns of an array property.

    The outputter writes a `ColumnValues` dataset only for the first property on a given grid of columns, so later
    properties on the same grid have none of their own. In that case the grid of the `fallback` property is used, which
    must be on the same grid.
    """
    if name + "ColumnValues" in nodes:
        return nodes[name + "ColumnValues"][:]
    return nodes[fallback + "ColumnValues"][:]


def heatingGrid(nodes, length):
    """Return the wavelength grid of the heating spectra, which all share one grid in this test."""
    grids = [nodes[name][:] for name in nodes.keys() if "StellarSED" in name and name.endswith("ColumnValues")]
    grids = [grid for grid in grids if grid.size == length]
    if not grids or any(not np.array_equal(grid, grids[0]) for grid in grids):
        print("FAILED: the heating spectra do not share a single wavelength grid, which this test assumes")
        sys.exit(0)
    return grids[0]


def absorbedLuminosities(nodes, attenuator):
    """Return the luminosity absorbed by each phase of the given attenuator, in L☉, for each galaxy."""
    absorbed = {}
    marker   = f":dustAbsorbed:{attenuator}:"
    for name in nodes.keys():
        if marker not in name or name.endswith("ColumnValues") or name.endswith("Columns"):
            continue
        phase = name.split(":")[-1]
        data  = nodes[name][:]
        # Units are recorded in a compound `units` attribute, whose `unitsInSI` field gives the conversion to SI.
        if "units" not in nodes[name].attrs:
            print(f"FAILED: no units attribute on '{name}', so its luminosity can not be converted")
            sys.exit(0)
        units = float(nodes[name].attrs["units"]["unitsInSI"])
        if data.ndim == 2:
            # An absorbed spectrum, L_ν: integrate over frequency by the trapezoidal rule.
            wavelengths = nodes[name + "ColumnValues"][:] if name + "ColumnValues" in nodes else heatingGrid(nodes, data.shape[1])
            frequencies = SPEED_OF_LIGHT / wavelengths
            luminosity  = np.sum(0.5 * (data[:, 1:] + data[:, :-1]) * np.abs(frequencies[:-1] - frequencies[1:]), axis=1)
        else:
            luminosity  = data
        absorbed[phase] = absorbed.get(phase, 0.0) + luminosity * units / LUMINOSITY_SOLAR
    return absorbed


def emittedLuminosity(nodes, name, fallback):
    """Return the integral over frequency of an emission spectrum, in L☉, summing over its resolution elements."""
    spectrum    = nodes[name][:]
    wavelengths = columnValues(nodes, name, fallback)
    return np.sum(spectrum * (SPEED_OF_LIGHT / wavelengths) * (FACTOR - 1.0 / FACTOR), axis=1)


with h5py.File(outputPath, "r") as f:
    # Outputs are numbered in order of time, so identify each by its expansion factor rather than assuming an order.
    expansionFactors = {name: float(f["Outputs"][name].attrs["outputExpansionFactor"]) for name in f["Outputs"].keys()}
    for outputName, expansionFactor in sorted(expansionFactors.items()):
        redshiftLabel = f"z={1.0 / expansionFactor - 1.0:.1f}"
        nodes = f["Outputs"][outputName]["nodeData"]
        for attenuator, phases in ATTENUATORS.items():
            label = f"{attenuator} at {redshiftLabel}"
            try:
                absorbed = absorbedLuminosities(nodes, attenuator)
                total    = f"dustEmissionSED:{attenuator}"
                emitted  = emittedLuminosity   (nodes, total, total)
                emittedPhases = {phase: emittedLuminosity(nodes, f"{total}:{phase}", total) for phase in phases}
            except KeyError as e:
                print(f"FAILED: {label}: expected dataset or attribute missing from the output: {e}")
                sys.exit(0)
            if sorted(absorbed.keys()) != sorted(phases):
                print(f"FAILED: {label}: absorbed luminosities found for phases {sorted(absorbed.keys())}, expected {sorted(phases)}")
                sys.exit(0)
            absorbedTotal = sum(absorbed.values())
            heated        = absorbedTotal > 0.0
            if not heated.any():
                print(f"FAILED: {label}: no galaxy absorbs any light, so energy balance is not being tested")
                sys.exit(0)
            worst = float(np.max(np.abs(emitted[heated] - absorbedTotal[heated]) / absorbedTotal[heated]))
            if worst > TOLERANCE_ENERGY:
                print(f"FAILED: {label}: emitted luminosity differs from absorbed luminosity by up to {worst:.3e}")
                sys.exit(0)
            print(f"SUCCESS: {label}: dust emits the luminosity it absorbs, to {worst:.3e}, for {int(heated.sum())} galaxies")
            for phase in phases:
                heatedPhase = absorbed[phase] > 0.0
                if not heatedPhase.any():
                    print(f"SUCCESS: {label}: the {phase} phase absorbs nothing, and emits {float(np.max(np.abs(emittedPhases[phase]))):.3e}")
                    continue
                worstPhase = float(np.max(np.abs(emittedPhases[phase][heatedPhase] - absorbed[phase][heatedPhase]) / absorbed[phase][heatedPhase]))
                if worstPhase > TOLERANCE_ENERGY:
                    print(f"FAILED: {label}: the {phase} phase emits a luminosity differing from that it absorbs by up to {worstPhase:.3e}")
                    sys.exit(0)
                print(f"SUCCESS: {label}: the {phase} phase emits the luminosity it absorbs, to {worstPhase:.3e}")

    # The observed frame. At redshift one, the observed-frame grid is exactly twice the rest-frame grid, so each bin is
    # evaluated at the same rest-frame wavelengths, and L_ν is not rescaled.
    outputRedshiftOne = [name for name, expansionFactor in expansionFactors.items() if abs(expansionFactor - 0.5) < 1.0e-6]
    if len(outputRedshiftOne) != 1:
        print(f"FAILED: expected exactly one output at redshift one, found expansion factors {sorted(expansionFactors.values())}")
        sys.exit(0)
    nodes = f["Outputs"][outputRedshiftOne[0]]["nodeData"]
    try:
        rest     = nodes["dustEmissionSED:charlotFall2000"         ][:]
        observed = nodes["dustEmissionSED:charlotFall2000:observed"][:]
        wavelengthsRest     = columnValues(nodes, "dustEmissionSED:charlotFall2000"         , "dustEmissionSED:charlotFall2000")
        wavelengthsObserved = columnValues(nodes, "dustEmissionSED:charlotFall2000:observed", "dustEmissionSED:charlotFall2000:observed")
    except KeyError as e:
        print(f"FAILED: expected observed-frame dataset missing from the output: {e}")
        sys.exit(0)
    if not np.allclose(wavelengthsObserved, 2.0 * wavelengthsRest, rtol=1.0e-12, atol=0.0):
        print("FAILED: the observed-frame wavelength grid is not twice the rest-frame grid")
        sys.exit(0)
    emitting = rest > 0.0
    if not emitting.any():
        print("FAILED: no galaxy emits at redshift one, so the observed frame is not being tested")
        sys.exit(0)
    worstFrame = float(np.max(np.abs(observed[emitting] - rest[emitting]) / rest[emitting]))
    if worstFrame > TOLERANCE_FRAME:
        print(f"FAILED: the observed-frame emission at redshift one differs from the rest-frame emission by up to {worstFrame:.3e}")
        sys.exit(0)
    print(f"SUCCESS: the observed-frame emission at redshift one equals the rest-frame emission, to {worstFrame:.3e}")

sys.exit(0)
