#!/usr/bin/env python3
"""Exercise legacy and hybrid bar prescriptions on the existing isolated tree.

Written with Codex assistance. Compiled-model validation is required before
this test, or the new prescription, is considered verified.
"""
import argparse
import os
from pathlib import Path
import subprocess
from xml.etree import ElementTree as ET

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def parameters(name, directory):
    tree = ET.parse(ROOT / "testSuite/parameters/barInstability.xml")
    root = tree.getroot()
    root.find("galacticStructureSolver").set("value", "equilibrium")
    root.find("outputFileName").set("value", str(directory / (name + ".hdf5")))
    bar = root.find("nodeOperator/nodeOperator[@value='barInstability']")
    if name != "legacy":
        model = ET.SubElement(bar, "barInstabilitySpheroidAngularMomentum",
                              value="retained" if name == "retained" else "bindingEnergy")
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
    names = ("legacy", "retained", "hybrid_default", "hybrid_explicit")
    paths = {name: parameters(name, directory) for name in names}
    if args.prepare_only:
        print("Prepared test inputs only; no compiled-model tests have run.")
        return
    fields = ("diskMassStellar", "diskMassGas", "spheroidMassStellar",
              "spheroidMassGas", "diskRadius", "spheroidRadius")
    results = {}
    for name, path in paths.items():
        with (directory / (name + ".log")).open("w") as log:
            subprocess.run([str(args.executable.resolve()), str(path)], cwd=ROOT,
                           env={**os.environ, "OMP_NUM_THREADS": "1"},
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        with h5py.File(directory / (name + ".hdf5")) as output:
            assert output.attrs["statusCompletion"] == 0, name
            data = output["Outputs/Output1/nodeData"]
            results[name] = {field: data[field][:] for field in fields}
        for values in results[name].values():
            assert np.all(np.isfinite(values)) and np.all(values >= 0), name
        spheroid = results[name]["spheroidMassStellar"] + results[name]["spheroidMassGas"]
        disk = results[name]["diskMassStellar"] + results[name]["diskMassGas"]
        assert np.any(spheroid > 0.01 * (spheroid + disk)), "Test did not build an established spheroid"
        assert np.all(results[name]["spheroidRadius"][spheroid > 0] > 0), name
    for left, right in (("legacy", "retained"), ("hybrid_default", "hybrid_explicit")):
        for field in fields:
            np.testing.assert_allclose(results[left][field], results[right][field], rtol=1e-9, atol=0,
                                       err_msg=f"{left} vs {right}: {field}")
    print("SUCCESS: bar prescription defaults, explicit selection, and finite spheroid properties")


if __name__ == "__main__":
    main()
