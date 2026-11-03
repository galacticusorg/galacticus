#!/usr/bin/env python3
"""Exercise the bar-instability spheroid angular-momentum prescriptions on the existing isolated tree.

Checks that:

* the default prescription and an explicit ``angularMomentumTransfer`` selection give identical results;
* the binding-energy prescription's documented transition defaults are its actual defaults;
* every variant builds an established spheroid with finite, non-negative properties;
* with converged numerical settings, the binding-energy prescription conserves its binding-energy proxy once the spheroid is
  established. In this fixture star formation only converts gas to stars within each component and there are no outflows, so,
  apart from a slow drift as the halo evolves, the bar is the only process that changes the proxy. As a control, the
  angular-momentum-transfer prescription, which does not conserve the proxy, must change it by much more than the tolerance
  over the same interval.
"""
import argparse
import os
from pathlib import Path
import subprocess
import sys
from xml.etree import ElementTree as ET

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

# Output times (Gyr) spanning the evolution of the isolated tree, which starts at 12.47 Gyr. Logarithmically spaced early times
# resolve the growth of the spheroid through the transition between prescriptions.
TIMES = np.concatenate((12.47 + np.geomspace(1.0e-5, 4.0e-3, 20), np.linspace(12.475, 13.47, 41)))
# Converged numerical settings for the conservation test, and the tolerance on the fractional change in the proxy.
STEP_CONVERGED = "3.1622776601683794e-4"
TOLERANCE_CONVERGED = "1.0e-7"
TOLERANCE_PROXY = 0.02
# The proxy's coefficients (the defaults) and the baryonic spheroid fraction above which the energy prescription applies fully.
FORM_FACTOR_DISK = 0.5
FORM_FACTOR_SPHEROID = 0.5
FORM_FACTOR_INTERACTION = 1.0
FRACTION_ESTABLISHED = 0.010
# Ratios of half-mass radius to scale radius for the exponential disk and the Hernquist spheroid.
HALF_MASS_DISK = 1.678346990016661
HALF_MASS_SPHEROID = 1.0 + np.sqrt(2.0)


def fail(message):
    print(f"FAILED: {message}")
    sys.exit(0)


def parameters(name, directory):
    tree = ET.parse(ROOT / "testSuite/parameters/barInstability.xml")
    root = tree.getroot()
    solver = root.find("galacticStructureSolver")
    solver.set("value", "equilibrium")
    root.find("outputFileName").set("value", str(directory / (name + ".hdf5")))
    outputTimes = ET.SubElement(root, "outputTimes", value="list")
    ET.SubElement(outputTimes, "times", value=" ".join(f"{time:.9f}" for time in TIMES))
    bar = root.find("nodeOperator/nodeOperator[@value='barInstability']")
    # The simple-size fixture uses 0.7 for stars, which leaves this disk stable
    # with equilibrium sizes. Use the usual stellar threshold in every variant.
    dynamics = bar.find("galacticDynamicsBarInstability")
    dynamics.find("stabilityThresholdStellar").set("value", "1.1")
    if name.endswith("_converged"):
        ET.SubElement(solver, "solutionTolerance", value=TOLERANCE_CONVERGED)
    if name not in ("legacy", "legacy_converged"):
        selector = "angularMomentumTransfer" if name == "angular_momentum_transfer" else "bindingEnergy"
        model = ET.SubElement(bar, "barInstabilitySpheroidAngularMomentum", value=selector)
        if name == "hybrid_explicit":
            ET.SubElement(model, "fractionBaryonicSpheroidTransitionMinimum", value="0.003")
            ET.SubElement(model, "fractionBaryonicSpheroidTransitionMaximum", value="0.010")
        if name == "energy_converged":
            ET.SubElement(model, "stepFiniteDifferenceRelative", value=STEP_CONVERGED)
    path = directory / (name + ".xml")
    tree.write(path, encoding="UTF-8", xml_declaration=True)
    return path


def proxy(series):
    """The binding-energy proxy (in units with G=1) and the baryonic spheroid fraction at each output."""
    massDisk = series["diskMassStellar"] + series["diskMassGas"]
    massSpheroid = series["spheroidMassStellar"] + series["spheroidMassGas"]
    radiusDisk = HALF_MASS_DISK * series["diskRadius"]
    radiusSpheroid = HALF_MASS_SPHEROID * series["spheroidRadius"]
    with np.errstate(divide="ignore", invalid="ignore"):
        energy = (+FORM_FACTOR_DISK * massDisk**2 / radiusDisk
                  + FORM_FACTOR_SPHEROID * massSpheroid**2 / radiusSpheroid
                  + FORM_FACTOR_INTERACTION * massDisk * massSpheroid / (radiusDisk + radiusSpheroid))
        fraction = massSpheroid / (massDisk + massSpheroid)
    return energy, fraction


def changeProxy(name, history):
    """Fractional change in the proxy from the first output with an established spheroid to the final output."""
    energy, fraction = history
    established = np.flatnonzero(fraction > FRACTION_ESTABLISHED)
    if len(established) < 10:
        fail(f"{name}: only {len(established)} outputs have an established spheroid; the conservation test needs at least 10")
    start = established[0]
    if not np.all(np.isfinite(energy[start:])):
        fail(f"{name}: binding-energy proxy is not finite after the spheroid is established")
    return energy[-1] / energy[start] - 1.0, TIMES[start]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, default=ROOT / "Galacticus.exe")
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    directory = ROOT / "testSuite/outputs/barBindingEnergy"
    directory.mkdir(parents=True, exist_ok=True)
    names = ("legacy", "angular_momentum_transfer", "hybrid_default", "hybrid_explicit", "legacy_converged", "energy_converged")
    paths = {name: parameters(name, directory) for name in names}
    if args.prepare_only:
        print("SKIPPED: prepared test inputs only; no compiled-model tests have run.")
        sys.exit(0)
    fields = ("diskMassStellar", "diskMassGas", "spheroidMassStellar",
              "spheroidMassGas", "diskRadius", "spheroidRadius")
    results = {}
    histories = {}
    for name, path in paths.items():
        print(f"Running {name}; log: {directory / (name + '.log')}", flush=True)
        with (directory / (name + ".log")).open("w") as log:
            status = subprocess.run([str(args.executable.resolve()), str(path)], cwd=ROOT,
                                    env={**os.environ, "OMP_NUM_THREADS": "1"},
                                    stdout=log, stderr=subprocess.STDOUT)
        if status.returncode != 0:
            fail(f"{name}: model failed to run (exit status {status.returncode}); see {directory / (name + '.log')}")
        try:
            with h5py.File(directory / (name + ".hdf5"), "r") as output:
                if output.attrs["statusCompletion"] != 0:
                    fail(f"{name}: model did not complete")
                outputs = sorted(output["Outputs"].keys(), key=lambda key: int(key[6:]))
                if len(outputs) != len(TIMES):
                    fail(f"{name}: expected {len(TIMES)} outputs, found {len(outputs)}")
                series = {field: np.array([output[f"Outputs/{key}/nodeData/{field}"][0] for key in outputs]) for field in fields}
        except (OSError, KeyError, IndexError) as error:
            fail(f"{name}: unable to read model output ({error})")
        for field, values in series.items():
            if not (np.all(np.isfinite(values)) and np.all(values >= 0)):
                fail(f"{name}: {field} is not finite and non-negative")
        results[name] = {field: values[-1:] for field, values in series.items()}
        histories[name] = proxy(series)
        spheroid = results[name]["spheroidMassStellar"] + results[name]["spheroidMassGas"]
        disk = results[name]["diskMassStellar"] + results[name]["diskMassGas"]
        fraction = np.divide(spheroid, spheroid + disk, out=np.zeros_like(spheroid), where=(spheroid + disk) > 0)
        print(f"{name}: final baryonic B/T = {fraction}", flush=True)
        if not np.any(spheroid > 0.01 * (spheroid + disk)):
            fail(f"{name}: test did not build an established spheroid; B/T={fraction}")
        if not np.all(results[name]["spheroidRadius"][spheroid > 0] > 0):
            fail(f"{name}: non-positive spheroid radius")
    for left, right in (("legacy", "angular_momentum_transfer"), ("hybrid_default", "hybrid_explicit")):
        for field in fields:
            if not np.allclose(results[left][field], results[right][field], rtol=1e-9, atol=0):
                fail(f"{left} vs {right}: {field} differs")
    change, timeStart = changeProxy("energy_converged", histories["energy_converged"])
    print(f"energy_converged: fractional change in the binding-energy proxy from t={timeStart:.6f} Gyr to the end = {change:+.4f}")
    if abs(change) > TOLERANCE_PROXY:
        fail(f"binding-energy prescription changed its proxy by {change:+.4f} after the spheroid was established (tolerance {TOLERANCE_PROXY})")
    changeControl, timeStart = changeProxy("legacy_converged", histories["legacy_converged"])
    print(f"legacy_converged: fractional change in the binding-energy proxy from t={timeStart:.6f} Gyr to the end = {changeControl:+.4f}")
    if abs(changeControl) < 5.0 * TOLERANCE_PROXY:
        fail(f"control: angular-momentum transfer changed the proxy by only {changeControl:+.4f}, so the conservation test cannot discriminate")
    print("SUCCESS: bar prescription defaults, angularMomentumTransfer selection, finite spheroid properties, and binding-energy conservation")


if __name__ == "__main__":
    main()
