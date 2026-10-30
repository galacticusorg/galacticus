#!/usr/bin/env python3
"""Exercise legacy and hybrid bar prescriptions on the existing isolated tree.

Written with Codex assistance. Compiled-model validation is required before
this test, or the new prescription, is considered verified.
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


def fail(message):
    print(f"FAILED: {message}")
    sys.exit(0)


def parameters(name, directory):
    tree = ET.parse(ROOT / "testSuite/parameters/barInstability.xml")
    root = tree.getroot()
    root.find("galacticStructureSolver").set("value", "equilibrium")
    root.find("outputFileName").set("value", str(directory / (name + ".hdf5")))
    bar = root.find("nodeOperator/nodeOperator[@value='barInstability']")
    # The simple-size fixture uses 0.7 for stars, which leaves this disk stable
    # with equilibrium sizes. Use the usual stellar threshold in every variant.
    dynamics = bar.find("galacticDynamicsBarInstability")
    dynamics.find("stabilityThresholdStellar").set("value", "1.1")
    if name != "legacy":
        selector = {"retained": "retained", "angular_momentum_transfer": "angularMomentumTransfer"}.get(name, "bindingEnergy")
        model = ET.SubElement(bar, "barInstabilitySpheroidAngularMomentum", value=selector)
        if name == "hybrid_explicit":
            ET.SubElement(model, "spheroidBaryonFractionTransitionMinimum", value="0.003")
            ET.SubElement(model, "spheroidBaryonFractionTransitionMaximum", value="0.010")
    path = directory / (name + ".xml")
    tree.write(path, encoding="UTF-8", xml_declaration=True)
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, default=ROOT / "Galacticus.exe")
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    directory = ROOT / "testSuite/outputs/barBindingEnergy"
    directory.mkdir(parents=True, exist_ok=True)
    names = ("legacy", "angular_momentum_transfer", "retained", "hybrid_default", "hybrid_explicit")
    paths = {name: parameters(name, directory) for name in names}
    if args.prepare_only:
        print("SKIPPED: prepared test inputs only; no compiled-model tests have run.")
        sys.exit(0)
    fields = ("diskMassStellar", "diskMassGas", "spheroidMassStellar",
              "spheroidMassGas", "diskRadius", "spheroidRadius")
    results = {}
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
                data = output["Outputs/Output1/nodeData"]
                results[name] = {field: data[field][:] for field in fields}
        except (OSError, KeyError) as error:
            fail(f"{name}: unable to read model output ({error})")
        for field, values in results[name].items():
            if not (np.all(np.isfinite(values)) and np.all(values >= 0)):
                fail(f"{name}: {field} is not finite and non-negative")
        spheroid = results[name]["spheroidMassStellar"] + results[name]["spheroidMassGas"]
        disk = results[name]["diskMassStellar"] + results[name]["diskMassGas"]
        fraction = np.divide(spheroid, spheroid + disk, out=np.zeros_like(spheroid), where=(spheroid + disk) > 0)
        print(f"{name}: baryonic B/T = {fraction}", flush=True)
        if not np.any(spheroid > 0.01 * (spheroid + disk)):
            fail(f"{name}: test did not build an established spheroid; B/T={fraction}")
        if not np.all(results[name]["spheroidRadius"][spheroid > 0] > 0):
            fail(f"{name}: non-positive spheroid radius")
    for left, right in (("legacy", "angular_momentum_transfer"), ("legacy", "retained"),
                        ("hybrid_default", "hybrid_explicit")):
        for field in fields:
            if not np.allclose(results[left][field], results[right][field], rtol=1e-9, atol=0):
                fail(f"{left} vs {right}: {field} differs")
    print("SUCCESS: bar prescription defaults, angularMomentumTransfer/retained compatibility, and finite spheroid properties")


if __name__ == "__main__":
    main()
